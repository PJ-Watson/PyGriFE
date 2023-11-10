import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSlider, QLabel, QComboBox, QColorDialog, QFrame, QCheckBox, QFileDialog, QPushButton
from PyQt6.QtGui import QPalette, QColor, QImage, QPixmap, QPainter
from QtImageViewer import QtImageViewer
from pathlib import Path
import astropy.io.fits as pf
from astropy.visualization import (LogStretch, AsinhStretch, ManualInterval, SqrtStretch, LinearStretch)
import qimage2ndarray
import numpy as np
import time

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("Object Selection")

        layout_h = QHBoxLayout()

        layout_side = QVBoxLayout()

        self.seg_dir = None
        dir_sel_button = QPushButton('Select Directory', self)
        dir_sel_button.clicked.connect(self.change_directory)
        layout_side.addWidget(dir_sel_button)

        # Image stretch
        self.stretch = SqrtStretch()
        stretch_box = QComboBox()
        stretch_box.addItems(["Logarithmic", "Square Root", "Asinh", "Linear"])
        stretch_box.setCurrentText("Square Root")
        stretch_box.currentTextChanged.connect( self.stretch_update )

        stretch_label = QLabel("Stretch:", self)
        layout_side.addWidget(stretch_label)
        layout_side.addWidget(stretch_box)

        # Interval
        self.interval = ManualInterval(0,1)
        interval_box = QComboBox()
        self.interval_keys = ["minmax", "99.9%", "99.8%", "99.5%", "99%", "98%", "95%"]
        interval_box.addItems(self.interval_keys)
        interval_box.setCurrentText("99.8%")
        interval_box.currentTextChanged.connect( self.interval_update )

        interval_label = QLabel("Interval:", self)
        layout_side.addWidget(interval_label)
        layout_side.addWidget(interval_box)

        # Opacity
        self.opacity = 0
        self.opacity_box = QComboBox()
        self.opacity_box.addItems(["100%", "90%", "75%", "50%", "25%", "0%"])
        self.opacity_box.setCurrentText("0%")
        self.opacity_box.currentTextChanged.connect( self.opacity_update )

        opacity_label = QLabel("Opacity:", self)
        layout_side.addWidget(opacity_label)
        layout_side.addWidget(self.opacity_box)

        self.seg_text = QLabel(f"x={0: <6}y={0: <6}\nID={None}", self)
        layout_side.addWidget(self.seg_text)

        self.selected_ids = []

        self.invert_box = QCheckBox("Invert image")
        self.invert_box.stateChanged.connect(self.opacity_update)
        self.invert_box.setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px;}")
        layout_side.addWidget(self.invert_box)

        layout_side.addWidget(QLabel("Background colour:", self))

        self.frm = QFrame(self)
        self.frm.setStyleSheet("QWidget { background-color: %s }"
                               % QColor(120, 120, 120).name())
        self.frm.mousePressEvent = self.choose_background_colour
        layout_side.addWidget(self.frm)

        layout_h.addLayout(layout_side)

        self.viewer = QtImageViewer()
        self.viewer.aspectRatioMode = Qt.AspectRatioMode.KeepAspectRatio
        self.viewer.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.viewer.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.viewer.regionZoomButton = Qt.MouseButton.LeftButton  # set to None to disable
        self.viewer.zoomOutButton = Qt.MouseButton.RightButton  # set to None to disable
        self.viewer.wheelZoomFactor = 1.1  # Set to None or 1 to disable

        self.viewer.setMouseTracking(True)
        self.viewer.mousePositionOnImageChanged.connect(self.seg_text_update)
        self.viewer.leftMouseButtonReleased.connect(self.click_location)

        self.current_seg_id = 0

        self.viewer.panButton = Qt.MouseButton.MiddleButton  
        layout_h.addWidget(self.viewer)

        widget = QWidget()
        widget.setLayout(layout_h)
        self.setCentralWidget(widget)
        self.change_directory()

    def change_directory(self, event=None):

        if self.seg_dir is None:
            self.seg_dir = Path.home()
        dir_name = QFileDialog.getExistingDirectory(self, "Open directory", str(self.seg_dir))
        if dir_name:
            self.seg_dir = Path(dir_name)
            try:
                self.load_image()
            except Exception as e:
                print (e)

    def click_location(self, modifiers, x, y):

        if not hasattr(self, "seg_data"):
            return
        x = int(x)
        y = int(y)
        seg_id = self.seg_data[y, x]

        if (len(self.selected_ids) > 0) & (modifiers == Qt.KeyboardModifier.ControlModifier):
            if seg_id == 0:
                pass
            elif self.selected_ids == [0]:
                self.selected_ids = [seg_id]
            else:
                self.selected_ids.append(seg_id)
        else:
            self.selected_ids = [seg_id]

        self.highlight_section(self.selected_ids)

    def highlight_section(self, seg_id):

        if seg_id == [0]:
            seg_map = np.zeros_like(self.seg_data, dtype="uint8")
        else:
            seg_map = np.isin(self.seg_data, seg_id).astype('uint8')*255
        seg_plot = self.seg_q.copy()
        seg_plot.setAlphaChannel(QImage(seg_map, seg_map.shape[1], seg_map.shape[0], seg_map.strides[0], QImage.Format.Format_Indexed8))
        pixmap = QPixmap.fromImage(seg_plot)

        if not hasattr(self.viewer, "overlay"):
            self.viewer.overlay = self.viewer.scene.addPixmap(pixmap)
        else:
            self.viewer.overlay.setPixmap(pixmap)

    def seg_text_update(self, pos):        
        
        if not hasattr(self, "seg_data"):
            return

        x = int(pos.x())
        y = int(pos.y())
        seg_id = self.seg_data[y, x]

        self.seg_text.setText(f"x={x: <6}y={y: <6}\nID={seg_id}")


    def interval_update(self, value):
        self.interval = ManualInterval(self.interval_dict[value][0], self.interval_dict[value][1])
        self.reload_image()

    def stretch_update(self, value):
        match value:
            case "Linear":
                self.stretch = LinearStretch()
            case "Asinh":
                self.stretch = AsinhStretch()
            case "Square Root":
                self.stretch = SqrtStretch()
            case "Logarithmic":
                self.stretch = LogStretch()
        self.reload_image()

    def opacity_update(self, event=None):
        
        self.opacity = float(self.opacity_box.currentText().split("%")[0])/100

        if self.invert_box.checkState() == Qt.CheckState.Checked:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.seg_mask+self.opacity_mask, a_min=0, a_max=1)*255
        else:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.opacity_mask+self.seg_mask, a_min=0, a_max=1)*255

        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)


    def reload_image(self):

        self.img_array[:,:,:-1] = (self.stretch(self.interval(self.data_array[:,:,:-1]))*255).astype('uint8')
        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)

    def load_image(self):

        t1 = time.time()
        print ("Reading images...", end="\r")
        
        with pf.open(self.seg_dir / "nis-wfss-ir_seg.fits") as hdul_seg:
            self.seg_mask = (hdul_seg[0].data > 0).astype('uint8')[::-1, :]
            self.opacity_mask = 1 - self.seg_mask
            self.seg_data = hdul_seg[0].data[::-1,:]
            self.seg_q = qimage2ndarray.array2qimage(np.stack(
                    [
                        np.zeros_like(self.seg_mask), 
                        np.ones_like(self.seg_mask), 
                        np.zeros_like(self.seg_mask), 
                        self.seg_mask*0.5,
                    ],
                    axis=-1), True)

        self.data_array = np.zeros((self.seg_mask.shape[0], self.seg_mask.shape[1], 4))
        self.overlap_mask = np.zeros_like(self.seg_mask, dtype="bool")

        self.data_array[:,:,-1] = self.seg_mask

        for i, filt in enumerate(["f200w", "f150w", "f115w"]):
            
            with pf.open(self.seg_dir / f"nis-wfss-{filt}_drz_sci.fits") as hdul_sci:
                self.data_array[:,:,i] = hdul_sci[0].data[::-1, :]
            with pf.open(self.seg_dir / f"nis-wfss-{filt}_drz_wht.fits") as hdul_wht:
                self.overlap_mask = self.overlap_mask | (hdul_wht[0].data[::-1, :] > 0)

        self.opacity_mask[~self.overlap_mask] = 0
        img_array = self.data_array.copy()
        print ("Reading images... DONE")

        print ("Computing intervals...", end="\r")
        self.interval_dict = {}
        percentiles = []
        for interval in self.interval_keys:
            match interval:
                case "minmax":
                    percentiles.append(100)
                case _:
                    percentiles.append(float(interval.split("%")[0]))
        lims = self.calc_interval_limits(percentiles)
        for k, l in zip(self.interval_keys, lims):
            self.interval_dict[k] = l
        
        self.interval = ManualInterval(self.interval_dict["99.8%"][0], self.interval_dict["99.8%"][1])
        print ("Computing intervals... DONE")

        print ("Formatting image for display...", end="\r")
        img_array[:,:,:-1] = self.stretch(self.interval(img_array[:,:,:-1]))
        self.img_array = (img_array*255).astype('uint8')
        self.q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(self.q_img)
        self.viewer.setBackgroundBrush(QColor(120,120,120))

        print ("Formatting image for display... DONE")
        print (f"Completed in {(time.time()-t1):.2f}s")

    def calc_interval_limits(self, percentiles):
        
        all_p = []
        for p in percentiles:
            lower_percent = (100 - p)*0.5
            upper_percent = 100 - lower_percent
            all_p.extend([lower_percent, upper_percent])
        limits = np.percentile(self.data_array[:,:,:-1].ravel(), all_p)

        res = [*zip(limits[::2], limits[1::2])]
        return res

    def choose_background_colour(self, event=None):

        col = QColorDialog.getColor()

        if col.isValid():

            self.frm.setStyleSheet("QWidget { background-color: %s }" 
                                   % col.name())
            self.viewer.setBackgroundBrush(col)


if __name__=="__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.showMaximized()
    app.exec()