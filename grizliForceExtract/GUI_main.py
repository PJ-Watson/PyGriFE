import sys
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSlider, QLabel, QComboBox, QColorDialog, QFrame, QCheckBox, QFileDialog, QPushButton, QFormLayout, QLineEdit
from PyQt6.QtGui import QPalette, QColor, QImage, QPixmap, QPainter
from QtImageViewer import QtImageViewer
from pathlib import Path
import astropy.io.fits as pf
from astropy.visualization import (LogStretch, AsinhStretch, ManualInterval, SqrtStretch, LinearStretch)
from astropy import wcs
import qimage2ndarray
import numpy as np
import time

class MainWindow(QMainWindow):
    def __init__(self, filters=["F115W", "F150W", "F200W"]):
        super(MainWindow, self).__init__()

        self.filters = filters

        self.setWindowTitle("Object Selection")

        layout_h = QHBoxLayout()

        left_toolbar = QWidget(self)
        left_toolbar.setFixedWidth(175)
        layout_side = QVBoxLayout(left_toolbar)
        layout_side.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.seg_dir = None
        self.files_window = None
        dir_sel_button = QPushButton('Select Files', self)
        dir_sel_button.clicked.connect(self.open_files_window)
        layout_side.addWidget(dir_sel_button)
        layout_side.addWidget(Separator())

        self.seg_text = QLabel(f"\n\n", self)
        layout_side.addWidget(self.seg_text)
        layout_side.addWidget(Separator())

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


        self.selected_ids = []

        self.invert_box = QCheckBox("Invert image")
        self.invert_box.stateChanged.connect(self.opacity_update)
        p = QPalette(self.invert_box.palette())
        p.setColor(
            QPalette.ColorGroup.Active, QPalette.ColorRole.Base, QColor(90,90,90),
            )
        self.invert_box.setPalette(p)
        layout_side.addWidget(self.invert_box)

        layout_side.addWidget(QLabel("Background colour:", self))

        self.bkg_frm = QFrame(self)
        self.bkg_frm.mousePressEvent = self.choose_background_colour
        self.bkg_frm.setMinimumHeight(50)
        self.bkg_frm.bkg_col = QColor("#787878")
        self.bkg_frm.setStyleSheet(
            f"QWidget {{ background-color: {self.bkg_frm.bkg_col.name()} }}"
        )
        layout_side.addWidget(self.bkg_frm)


        combine_button = QPushButton('Combine Selection', self)
        combine_button.clicked.connect(self.combine_ids)
        save_button = QPushButton('Save Map', self)
        save_button.clicked.connect(self.save_output)
        layout_side.addWidget(Separator())
        layout_side.addWidget(combine_button)
        layout_side.addWidget(save_button)

        layout_h.addWidget(left_toolbar)

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

        # self.change_directory()

    def open_files_window(self, event=None):

        if self.files_window is None:
            self.files_window = FilesWindow(self)
        self.files_window.show()

    def change_directory(self, event=None):

        if self.seg_dir is None:
            self.seg_dir = Path.home()
        dir_name = QFileDialog.getExistingDirectory(self, "Open directory", str(self.seg_dir))
        if dir_name:
            self.seg_dir = Path(dir_name)
            # try:
            #     self.load_image()
            # except Exception as e:
            #     print (e)

    def click_location(self, modifiers, x, y):

        if not hasattr(self, "img_array"):
            return
        x = int(x)
        y = int(y)
        seg_id = self.seg_data[y, x]

        if seg_id in self.selected_ids:
            self.selected_ids.remove(seg_id)
        elif (len(self.selected_ids) > 0) & (modifiers == Qt.KeyboardModifier.ControlModifier):
            if seg_id == 0:
                pass
            elif self.selected_ids == [0]:
                self.selected_ids = [seg_id]
            else:
                self.selected_ids.append(seg_id)
        else:
            self.selected_ids = [seg_id]

        if self.selected_ids == []:
            self.selected_ids = [0]

        self.highlight_section(self.selected_ids)

    def highlight_section(self, seg_id):

        if not hasattr(self, "img_array"):
            return

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
        
        if not hasattr(self, "img_array"):
            return

        x = int(pos.x())
        y = int(pos.y())
        seg_id = self.seg_data[y, x]

        self.seg_text.setText(
            f"{wcs.utils.pixel_to_skycoord(x, y, self.wcs).to_string(precision=6)}\n"
            f"x={x: <6}y={y: <6}\nID={seg_id}"
        )


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

        if not hasattr(self, "img_array"):
            return

        if self.invert_box.checkState() == Qt.CheckState.Checked:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.seg_mask+self.opacity_mask, a_min=0, a_max=1)*255
        else:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.opacity_mask+self.seg_mask, a_min=0, a_max=1)*255

        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)


    def reload_image(self):

        if not hasattr(self, "img_array"):
            return

        self.img_array[:,:,:-1] = (self.stretch(self.interval(self.data_array[:,:,:-1]))*255).astype('uint8')
        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)

    def load_image(self):

        t1 = time.time()
        print ("Reading images...", end="\r")
        
        with pf.open(self.seg_img_path) as hdul_seg:
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

            self.wcs = wcs.WCS(hdul_seg[0].header)

        self.data_array = np.zeros((self.seg_mask.shape[0], self.seg_mask.shape[1], 4))
        self.overlap_mask = np.zeros_like(self.seg_mask, dtype="bool")

        self.data_array[:,:,-1] = self.seg_mask

        for i, filt_img_path in enumerate([self.r_img_path, self.g_img_path, self.b_img_path]):
            try:
                with pf.open(filt_img_path) as hdul_sci:
                    self.data_array[:,:,i] = hdul_sci[0].data[::-1, :]
                try:
                    wht_path = filt_img_path.replace("sci", "wht")
                    with pf.open(wht_path) as hdul_wht:
                        self.overlap_mask = self.overlap_mask | (hdul_wht[0].data[::-1, :] > 0)
                except:
                    pass
            except:
                pass

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

        self.opacity = float(self.opacity_box.currentText().split("%")[0])/100

        if self.invert_box.checkState() == Qt.CheckState.Checked:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.seg_mask+self.opacity_mask, a_min=0, a_max=1)*255
        else:
            self.img_array[:,:,-1] = np.clip(self.opacity*self.opacity_mask+self.seg_mask, a_min=0, a_max=1)*255

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

        col = QColorDialog.getColor(self.bkg_frm.bkg_col)

        if col.isValid():
                
            self.bkg_frm.bkg_col = col
            self.bkg_frm.setStyleSheet(
                f"QWidget {{ background-color: {self.bkg_frm.bkg_col.name()} }}"
            )
            self.viewer.setBackgroundBrush(col)

    def combine_ids(self, event=None):
        if not hasattr(self, "seg_data"):
            print ("No segmentation mask to modify.")
            return

        current_selection = np.isin(self.seg_data, self.selected_ids)

        self.seg_data[current_selection] = np.nanmin(self.selected_ids)

    def save_output(self, event=None):

        if not (hasattr(self, "seg_img_path") and hasattr(self, "seg_data")):
            print ("No segmentation mask loaded.")
            return

        backup_path = self.seg_img_path.parent / f"{self.seg_img_path.stem}_backup.fits"
        with pf.open(self.seg_img_path) as seg_hdul:
            if not backup_path.is_file():
                seg_hdul.writeto(backup_path)
            
            seg_hdul[0].data = self.seg_data[::-1,:]

            seg_hdul.writeto(self.seg_img_path.parent / f"{self.seg_img_path.stem}_modified.fits")

class Separator(QFrame):
    def __init__(self):
        super(QFrame, self).__init__()
        self.setFrameShape(QFrame.Shape.HLine)
        self.setFrameShadow(QFrame.Shadow.Sunken)
        self.setLineWidth(3)

class cQLineEdit(QLineEdit):
    clicked=pyqtSignal()
    def __init__(self,*args, **kwargs):
        super().__init__(*args,**kwargs)

    def mousePressEvent(self,QMouseEvent):
        # print (self.text())
        # print (self.parent())
        if (self.text is None or self.text=="") and (self.parent().recent_dir is None):
            init = str(Path.home())
        else:
            init = self.text()
        f, _ = QFileDialog.getOpenFileName(self, "Select File", init, "FITS files (*.fits)")

        if f:
            self.setText(f)
            self.parent().recent_dir = Path(f).parent
            # print (self.parent().recent_dir)
            # try:
                # Path(f)
        # .getExistingDirectory(self, "Open directory", str(self.seg_dir))
        # self.clicked.emit()

class FilesWindow(QWidget):
    def __init__(self, root):
        super().__init__()
        self.root = root

        layout = QVBoxLayout()
        # self.label = QLabel("Another Window ")
        # layout.addWidget(self.label)

        self.recent_dir = None

        dir_sel_button = QPushButton('Fill From Directory', self)
        dir_sel_button.clicked.connect(self.change_directory)
        layout.addWidget(dir_sel_button)
        layout.addWidget(Separator())

        sub_layout = QFormLayout()
        self.seg_line = cQLineEdit(self)
        sub_layout.addRow("Segmentation Map:", self.seg_line)
        sub_layout.addRow(Separator())
        self.stack_line = cQLineEdit(self)
        sub_layout.addRow("Stacked Image:", self.stack_line)
        sub_layout.addRow(Separator())
        self.b_line = cQLineEdit(self)
        sub_layout.addRow("Blue:", self.b_line)
        self.g_line = cQLineEdit(self)
        sub_layout.addRow("Green:", self.g_line)
        self.r_line = cQLineEdit(self)
        sub_layout.addRow("Red:", self.r_line)
        layout.addLayout(sub_layout)

        load_all_button = QPushButton('Load Images', self)
        load_all_button.clicked.connect(self.load_all)
        layout.addWidget(Separator())
        layout.addWidget(load_all_button)

        self.setLayout(layout)
        self.setMinimumWidth(540)

    def change_directory(self, event=None):

        if self.root.seg_dir is not None:
            init = str(seg_dir)
        elif self.recent_dir is None:
            init = str(self.recent_dir)
        else:
            init = str(Path.home())

        dir_name = QFileDialog.getExistingDirectory(self, "Open directory", init)
        if dir_name:
            self.seg_dir = Path(dir_name)
            print (self.seg_dir)

            try:
                self.seg_line.setText(str([*self.seg_dir.glob("*ir_seg.fits")][0]))
            except:
                print ("Segmentation map not found.")

            try:
                self.stack_line.setText(str([*self.seg_dir.glob("*ir_drz_sci.fits")][0]))
            except:
                print ("Stacked image not found.")

            for f, l in zip(self.root.filters, [self.b_line, self.g_line, self.r_line]):
                l.setText(str([*self.seg_dir.glob(f"*{f.lower()}_drz_sci.fits")][0]))

    def load_all(self):
        self.root.seg_img_path = Path(self.seg_line.text())
        self.root.b_img_path = Path(self.b_line.text())
        self.root.g_img_path = Path(self.g_line.text())
        self.root.r_img_path = Path(self.r_line.text())
        self.root.load_image()
        

    def printText(self):
        print("This works")
    
if __name__=="__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.showMaximized()
    app.exec()

# import sys
# import time
# from PyQt6.QtWidgets import *
# from PyQt6.QtGui import *
# from PyQt6.QtCore import *


# class ProcessRunnable(QRunnable):
#     def __init__(self, target, args):
#         QRunnable.__init__(self)
#         self.t = target
#         self.args = args

#     def run(self):
#         self.t(*self.args)

#     def start(self):
#         QThreadPool.globalInstance().start(self)


# def run(user_input, log):
#     text = ""
#     if user_input == "":
#         text = "Please enter a value\n"
#     else:
#         text = "Test"
#         # Sleep for 5 seconds
#         time.sleep(5)

#     QMetaObject.invokeMethod(log, "append", Qt.ConnectionType.QueuedConnection, Q_ARG(str, text))


# class LogginOutput(QTextEdit):
#     def __init__(self, parent=None):
#         super(LogginOutput, self).__init__(parent)

#         self.setReadOnly(True)

#         self.setLineWrapMode(self.LineWrapMode.NoWrap)

#         self.insertPlainText("")

#     @pyqtSlot(str)
#     def append(self, text):
#         self.moveCursor(QTextCursor.MoveOperation.End)
#         current = self.toPlainText()

#         if current == "":
#             self.insertPlainText(text)
#         else:
#             self.insertPlainText("\n" + text)

#         sb = self.verticalScrollBar()
#         sb.setValue(sb.maximum())


# class App(QWidget):
#     def __init__(self):
#         super().__init__()

#         self.init_ui()

#     def init_ui(self):
#         label = QLabel("Amount")
#         amount_input = QLineEdit()
#         submit = QPushButton("Submit", self)
#         box = LogginOutput(self)

#         submit.clicked.connect(lambda: self.changeLabel(box, amount_input))

#         grid = QGridLayout()
#         grid.addWidget(label, 0, 0)
#         grid.addWidget(amount_input, 1, 0)
#         grid.addWidget(submit, 1, 1)
#         grid.addWidget(box, 2, 0, 5, 2)

#         self.setLayout(grid)
#         self.resize(350, 250)
#         self.setWindowTitle("GetMeStuff Bot v0.1")
#         self.show()

#     def center(self):
#         qr = self.frameGeometry()
#         cp = QGuiApplication.primaryScreen().availableGeometry().center()
#         qr.moveCenter(cp)
#         self.move(qr.topLeft())

#     def changeLabel(self, box, user_input):
#         self.p = ProcessRunnable(target=run, args=(user_input.displayText(), box))
#         self.p.start()
#         user_input.clear()


# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     widget = App()
#     app.exec()

#     # app = QApplication(sys.argv)
#     # window = MainWindow()
#     # window.showMaximized()
#     # app.exec()