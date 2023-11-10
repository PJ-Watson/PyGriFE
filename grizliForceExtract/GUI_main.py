# from PyQt6 import QtWidgets

# # Only needed for access to command line arguments
# import sys

# # # You need one (and only one) QApplication instance per application.
# # # Pass in sys.argv to allow command line arguments for your app.
# # # If you know you won't use command line arguments QApplication([]) works too.
# # app = QApplication(sys.argv)

# # # Create a Qt widget, which will be our window.
# # window = QWidget()
# # window.show()  # IMPORTANT!!!!! Windows are hidden by default.

# # # Start the event loop.
# # app.exec()


# # # Your application won't reach here until you exit and the event
# # # loop has stopped.

# class MainWindow(QtWidgets.QMainWindow):
#     def __init__(self):
#         QtWidgets.QMainWindow.__init__(self)


# if __name__=="__main__":
#     app = QtWidgets.QApplication(sys.argv)
#     window = MainWindow()
#     window.show()
#     app.exec()

import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSlider, QLabel, QComboBox, QColorDialog, QFrame
from PyQt6.QtGui import QPalette, QColor, QImage, QPixmap, QPainter
from QtImageViewer import QtImageViewer
from pathlib import Path
import astropy.io.fits as pf
from astropy.visualization import (LogStretch, AsinhStretch, ManualInterval, SqrtStretch, LinearStretch)
import qimage2ndarray
import numpy as np
import time


# Custom slot for handling mouse clicks in our viewer.
# This example just prints the (row, column) matrix index
# of the image pixel that was clicked on.
def handleLeftClick(x, y):
    row = int(y)
    column = int(x)
    print("Pixel (row="+str(row)+", column="+str(column)+")")

def viewSegMap(viewer, input):
    seg_dir = Path("/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/Prep")
    with pf.open(seg_dir / "nis-wfss-ir_seg.fits") as hdul:
        # img = hdul[0].data
        interval = MinMaxInterval()
        stretch = LogStretch()
        img = qimage2ndarray.array2qimage(stretch(interval(hdul[0].data)), normalize=True)
    viewer.setImage(img)

# class ColourFrame(QWidget):

#     def __init__(self, colour):
#         super(ColourFrame, self).__init__()
#         self.setAutoFillBackground(True)

#         palette = self.palette()
#         palette.setColor(QPalette.ColorRole.Window, QColor(colour))
#         self.setPalette(palette)

#     def set_colour(self, qcolour):
#         # self.setStyleSheet("QWidget { background-color: %s }" 
#         #                            % colour)
        
#         self.setPalette(self.palette().setColor(QPalette.ColorRole.Window, qcolour))

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("Object Selection")

        layout_h = QHBoxLayout()

        layout_side = QVBoxLayout()
        # self.test_col = ColourFrame("green")
        # layout_side.addWidget(self.test_col)
        # layout_side.addWidget(Color('blue'))

        # self.direct_interval = MinMaxInterval()

        # Image stretch
        self.stretch = AsinhStretch()
        stretch_box = QComboBox()
        stretch_box.addItems(["Logarithmic", "Square Root", "Asinh", "Linear"])
        stretch_box.setCurrentText("Asinh")
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
        opacity_box = QComboBox()
        opacity_box.addItems(["100%", "90%", "75%", "50%", "25%", "0%"])
        opacity_box.setCurrentText("0%")
        opacity_box.currentTextChanged.connect( self.opacity_update )

        opacity_label = QLabel("Opacity:", self)
        layout_side.addWidget(opacity_label)
        layout_side.addWidget(opacity_box)

        self.seg_text = QLabel(f"x={0: <6}y={0: <6}\nID={None}", self)
        layout_side.addWidget(self.seg_text)

        self.selected_ids = []

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

        
        # Allow panning with the middle mouse button.
        self.viewer.panButton = Qt.MouseButton.MiddleButton  # set to None to disable
            
        # Load an image file to be displayed (will popup a file dialog).
        # viewer.open()
        self.load_image()
        
        # viewer.leftMouseButtonReleased.connect(lambda input: viewSegMap(viewer, input))

        layout_h.addWidget(self.viewer)

        widget = QWidget()
        widget.setLayout(layout_h)
        self.setCentralWidget(widget)

    def click_location(self, modifiers, x, y):

        x = int(x)
        y = int(y)
        seg_id = self.seg_data[y, x]

        print (modifiers)

        if (len(self.selected_ids) > 0) & (modifiers == Qt.KeyboardModifier.ControlModifier):
            print (seg_id)
            if seg_id == 0:
                pass
            elif self.selected_ids == [0]:
                self.selected_ids = [seg_id]
            else:
                self.selected_ids.append(seg_id)
        else:
            self.selected_ids = [seg_id]

        print (self.selected_ids)

        self.highlight_section(self.selected_ids)

    def highlight_section(self, seg_id):

        import time
        t1 = time.time()

        print (seg_id)
        if seg_id == [0]:
            seg_map = np.zeros_like(self.seg_data, dtype="uint8")
        else:
            seg_map = np.isin(self.seg_data, seg_id).astype('uint8')*255
        # print (time.time() -t1)
        seg_plot = self.seg_q.copy()
        # QImage(seg_map, seg_map.shape[1], seg_map.shape[0], seg_map.strides[0], QImage.Format.Format_Indexed8)
        seg_plot.setAlphaChannel(QImage(seg_map, seg_map.shape[1], seg_map.shape[0], seg_map.strides[0], QImage.Format.Format_Indexed8))
        # seg_plot.setAlphaChannel(qimage2ndarray.gray2qimage(seg_map,True))
        pixmap = QPixmap.fromImage(seg_plot)

        if not hasattr(self.viewer, "overlay"):
            print ("Creating overlay")
            self.viewer.overlay = self.viewer.scene.addPixmap(pixmap)
        else:
            print ("modifying")
            self.viewer.overlay.setPixmap(pixmap)

        print (time.time() -t1)

            


    def seg_text_update(self, pos):

        x = int(pos.x())
        y = int(pos.y())
        seg_id = self.seg_data[y, x]

        self.seg_text.setText(f"x={x: <6}y={y: <6}\nID={seg_id}")


    def interval_update(self, value):
        match value:
            case "minmax":
                self.calc_interval_limits(100)
            case _:
                self.calc_interval_limits(float(value.split("%")[0]))
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
            # case "Sinh":
            #     self.stretch = SinhStretch()
            # case "Hist. Eq.":
            #     self.stretch = HistEqStretch()
            # case "Contrast Bias":
            #     self.stretch = ContrastBiasStretch()
                
        self.reload_image()

    def opacity_update(self, value):
        self.opacity = float(value.split("%")[0])/100

        self.img_array[:,:,-1] = np.clip(self.opacity*(1-self.seg_mask)+self.seg_mask, a_min=0, a_max=1)*255

        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)
        # self.reload_image()


    def reload_image(self):

        # img_array = np.stack(
        #     [
        #         self.stretch(self.interval(self.data_array[:,:,0])), 
        #         np.clip(self.opacity*(1-self.seg_mask)+self.seg_mask, a_min=0, a_max=1)
        #         #self.data_array[:,:,-1]
        #     ],
        #     axis=-1,
        # )
        # img_array = 
        t1 = time.time()
        # arr = self.data_array.copy()
        # print (img_array)
        print ("Copied:", time.time()- t1)
        # print (self.interval(self.data_array[:,:,:-1]))
        # self.calc_interval_limits()
        i1 = self.interval(self.data_array[:,:,:-1])
        print (time.time()-t1)
        i2 = self.stretch(i1)
        print (time.time()-t1)
        i2*=255
        print (time.time()-t1)
        self.img_array[:,:,:-1] = i2.astype('uint8')
        # self.img_array[:,:,:-1] = (self.stretch(self.interval(self.data_array[:,:,:-1]))*255).astype('uint8')
        print ("Stretched:", time.time()- t1)
        # print (self.interval.get_limits(img_array[:,:,:-1]))
        # arr[:,:,-1] = self.img_array[:,:,-1]
        # self.img_array[:,:,-1] = 
        # print ("Opacity:", time.time()- t1)
        # print (img_array[3000:3010,3000:3010])
        # img_array = (img_array*255).astype('uint8')
        # print (img_array[3000:3010,3000:3010])
        # img_array = (img_array*255).astype('uint8')
        print ("Multiplied:", time.time()- t1)
        q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        self.viewer.setImage(q_img)
        print ("Set image:", time.time()- t1)
        # self.alpha_channel = qimage2ndarray.alpha_view(self.viewer.image())
        # print (self.alpha_channel)

    def load_image(self):

        t1 = time.time()
        
        seg_dir = Path("/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/Prep")
        with pf.open(seg_dir / "nis-wfss-ir_seg.fits") as hdul_seg:
            self.seg_mask = (hdul_seg[0].data > 0).astype('uint8')[::-1, :]
            self.seg_data = hdul_seg[0].data[::-1,:]#.astype(int)
            # print (self.seg_data.dtype)
            # self.seg_data = self.seg_data [3000:4000,3000:4000]
            # self.seg_mask = self.seg_mask[3000:4000,3000:4000]

            # self.seg_overlay = np.stack(
            #         [
            #             np.zeros_like(self.seg_mask), 
            #             np.ones_like(self.seg_mask), 
            #             np.zeros_like(self.seg_mask), 
            #             self.seg_mask*0.5,
            #         ],
            #         axis=-1)
            self.seg_q = qimage2ndarray.array2qimage(np.stack(
                    [
                        np.zeros_like(self.seg_mask), 
                        np.ones_like(self.seg_mask), 
                        np.zeros_like(self.seg_mask), 
                        self.seg_mask*0.5,
                    ],
                    axis=-1), True)
            # print (np.sum(direct_data[3000:3100,3000:3100]))
            # print (np.sum(direct_data))
        # with pf.open(seg_dir / "nis-wfss-ir_drz_sci.fits") as hdul_drz:
        #     print ("here")
        #     direct_data = hdul_drz[0].data
        #     # direct_data = direct_data [3000:4000,3000:4000]
        #     # print (direct_data)
        #     self.data_array = np.stack([direct_data,direct_data,direct_data,self.seg_mask], axis=-1)

        self.data_array = np.zeros((self.seg_mask.shape[0], self.seg_mask.shape[1], 4))
        self.overlap_mask = np.ones_like(self.seg_mask, dtype="bool")

        self.data_array[:,:,-1] = self.seg_mask

        for i, filt in enumerate(["f200w", "f150w", "f115w"]):
            
            with pf.open(seg_dir / f"nis-wfss-{filt}_drz_sci.fits") as hdul_sci:
                self.data_array[:,:,i] = hdul_sci[0].data[::-1, :]
                print (np.percentile(hdul_sci[0].data[::-1, :].ravel(), (0.1, 99.9)))
            with pf.open(seg_dir / f"nis-wfss-{filt}_drz_wht.fits") as hdul_wht:
                self.overlap_mask = self.overlap_mask & (hdul_wht[0].data[::-1, :] > 0)

        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # plt.imshow(self.overlap_mask, alpha=0.5)
        # plt.imshow(np.log10(self.data_array[:,:,0]))
        # plt.show()
        # plt.imshow
        # print (np.sum(self.overlap_mask))
        # with pf.open(seg_dir / "nis-wfss-f150w_drz_sci.fits") as hdul_g:
        #     print ("here")
        #     g_data = hdul_g[0].data
        # with pf.open(seg_dir / "nis-wfss-f200w_drz_sci.fits") as hdul_r:
        #     print ("here")
        #     r_data = hdul_r[0].data
        #     # direct_data = direct_data [3000:4000,3000:4000]
        #     # print (direct_data)
        #     self.data_array = np.stack([r_data,g_data,b_data,self.seg_mask], axis=-1)

        
        # with pf.open(seg_dir / "nis-wfss-ir_drz_wht.fits") as hdul_wht:
        #     print ("wht")
        #     wht_data = hdul_wht[0].data
        #     # direct_data = direct_data [3000:4000,3000:4000]
        #     # print (direct_data)
        #     # print (self.data_array.shape)
        #     for i in range(3):
        #         # print (i.shape)
        #         self.data_array[:,:,i][wht_data==0] = np.nan


        img_array = self.data_array.copy()
        # print (img_array[:,:,:-1][self.overlap_mask].shape)
        # print (np.percentile(img_array[:,:,:-1][self.overlap_mask].ravel(), 99.95))

        self.calc_interval_limits(99.8)

        # print ("Computing intervals...")
        # for i in self.interval_keys:
        #     match i:
        #         case "minmax":
        #             self.calc_interval_limits(100)
        #         case _:
        #             self.calc_interval_limits(float(i.split("%")[0]))
        # print ("Computing intervals... DONE")
        img_array[:,:,:-1] = self.stretch(self.interval(img_array[:,:,:-1]))
        # print (img_array[3000:3010,3000:3010])
        self.img_array = (img_array*255).astype('uint8')
        # print (img_array[3000:3010,3000:3010])
        # print (img_array)

        # self.q_seg_img = qimage2ndarray.array2qimage(np.stack([np.full_like(self.seg_mask, 10000), self.seg_mask], axis=-1), normalize=True,)
        # self.q_seg_view = qimage2ndarray.recarray_view(self.q_seg_img)

        # self.q_img = qimage2ndarray.array2qimage(img_array, normalize=True)
        # print (f"Loaded in {time.time()-t1}")
        self.q_img = QImage(self.img_array, self.img_array.shape[1], self.img_array.shape[0], self.img_array.strides[0], QImage.Format.Format_RGBA8888)
        print (f"Loaded in {time.time()-t1}")
        self.viewer.setImage(self.q_img)

            # self.calc_interval_limits()

        self.viewer.setBackgroundBrush(QColor(120,120,120))
        # col = QColorDialog.getColor()
        # self.q_img_arr = qimage2ndarray.recarray_view(self.q_img)
        # self.alpha_channel = qimage2ndarray.alpha_view(self.viewer.image())
        # print (self.q_img_arr)

    def calc_interval_limits(self, percentile):

        lower_percent = (100 - percentile)*0.5
        upper_percent = 100 - lower_percent

        # print (self.data_array[self.overlap_mask].shape)
        # print (self.data_array[:,:,:-1][self.overlap_mask].shape)
        self.interval_vmin, self.interval_vmax = np.percentile(self.data_array[:,:,:-1][self.overlap_mask].ravel(), (lower_percent, upper_percent))
        # self.interval_vmin = -0.005
        # self.interval_vmax = 0.3
        self.interval = ManualInterval(vmin = self.interval_vmin, vmax = self.interval_vmax)
        # print (self.interval_vmin, self.interval_vmax)


    def choose_background_colour(self, event=None):

        col = QColorDialog.getColor()

        if col.isValid():

            self.frm.setStyleSheet("QWidget { background-color: %s }" 
                                   % col.name())
            self.viewer.setBackgroundBrush(col)


if __name__=="__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()

from PyQt6.QtWidgets import (QWidget, QPushButton, QFrame,
        QColorDialog, QApplication)
from PyQt6.QtGui import QColor
import sys


class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):

        col = QColor(0, 0, 0)

        self.btn = QPushButton('Dialog', self)
        self.btn.move(20, 20)

        self.btn.clicked.connect(self.showDialog)

        self.frm = QFrame(self)
        self.frm.setStyleSheet("QWidget { background-color: %s }"
                               % col.name())
        self.frm.setGeometry(130, 22, 200, 200)

        self.setGeometry(300, 300, 450, 350)
        self.setWindowTitle('Color dialog')
        self.show()


    def showDialog(self):

        col = QColorDialog.getColor()

        if col.isValid():

            self.frm.setStyleSheet("QWidget { background-color: %s }" 
                                   % col.name())


def main():

    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()


# if __name__ == '__main__':
#     # Create the QApplication.
#     app = QApplication(sys.argv)

#     seg_dir = Path("/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/Prep")

#     with pf.open(seg_dir / "nis-wfss-ir_drz_sci.fits") as hdul:
#         # img = hdul[0].data
#         interval = MinMaxInterval()
#         stretch = LogStretch()
#         img = qimage2ndarray.array2qimage(stretch(interval(hdul[0].data)), normalize=True)
        
#     # Create an image viewer widget.
#     viewer = QtImageViewer()
        
#     # Set viewer's aspect ratio mode.
#     # !!! ONLY applies to full image view.
#     # !!! Aspect ratio always ignored when zoomed.
#     #   Qt.AspectRatioMode.IgnoreAspectRatio: Fit to viewport.
#     #   Qt.AspectRatioMode.KeepAspectRatio: Fit in viewport using aspect ratio.
#     #   Qt.AspectRatioMode.KeepAspectRatioByExpanding: Fill viewport using aspect ratio.
#     viewer.aspectRatioMode = Qt.AspectRatioMode.KeepAspectRatio
    
#     # Set the viewer's scroll bar behaviour.
#     #   Qt.ScrollBarPolicy.ScrollBarAlwaysOff: Never show scroll bar.
#     #   Qt.ScrollBarPolicy.ScrollBarAlwaysOn: Always show scroll bar.
#     #   Qt.ScrollBarPolicy.ScrollBarAsNeeded: Show scroll bar only when zoomed.
#     viewer.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
#     viewer.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    
#     # Allow zooming by draggin a zoom box with the left mouse button.
#     # !!! This will still emit a leftMouseButtonReleased signal if no dragging occured,
#     #     so you can still handle left mouse button clicks in this way.
#     #     If you absolutely need to handle a left click upon press, then
#     #     either disable region zooming or set it to the middle or right button.
#     viewer.regionZoomButton = Qt.MouseButton.LeftButton  # set to None to disable
    
#     # Pop end of zoom stack (double click clears zoom stack).
#     viewer.zoomOutButton = Qt.MouseButton.RightButton  # set to None to disable
    
#     # Mouse wheel zooming.
#     viewer.wheelZoomFactor = 1.25  # Set to None or 1 to disable
    
#     # Allow panning with the middle mouse button.
#     viewer.panButton = Qt.MouseButton.MiddleButton  # set to None to disable
        
#     # Load an image file to be displayed (will popup a file dialog).
#     # viewer.open()
#     viewer.setImage(img)
    
#     # Handle left mouse clicks with your own custom slot
#     # handleLeftClick(x, y). (x, y) are image coordinates.
#     # For (row, col) matrix indexing, row=y and col=x.
#     # QtImageViewer also provides similar signals for
#     # left/right mouse button press, release and doubleclick.
#     # Here I bind the slot to leftMouseButtonReleased only because
#     # the leftMouseButtonPressed signal will not be emitted due to
#     # left clicks being handled by the regionZoomButton.
#     # viewer.leftMouseButtonReleased.connect(handleLeftClick)
#     viewer.leftMouseButtonReleased.connect(lambda input: viewSegMap(viewer, input))
        
#     # Show the viewer and run the application.
#     viewer.show()
#     sys.exit(app.exec())

# import sys
# from PyQt6.QtCore import Qt
# from PyQt6.QtWidgets import QApplication
# from QtImageStackViewer import QtImageStackViewer

# import numpy as np
    

# if __name__ == '__main__':

#     import matplotlib.pyplot as plt
#     arr = np.ones((5,5,3))

#     arr[0,1,0] = 0
#     arr[2,1,1] = 0
#     arr[4,2,2] = 0
#     arr[3,4,2] = 0

#     fig, axs = plt.subplots(1,4)
#     for i, a in enumerate(axs[:2]):
#         a.imshow(arr[:,:,i])

#     axs[3].imshow(arr[arr > 0].axis())
#     plt.show()

#     # Create the QApplication.
#     app = QApplication(sys.argv)

#     rng = np.random.default_rng()
#     arr_img = rng.random((5,5,3))
        
#     # Create an image stack viewer widget.
#     viewer = QtImageStackViewer(arr_img)
    
#     # Customize mouse interaction via the QtImageViewer widget.
#     viewer.imageViewer.regionZoomButton = Qt.MouseButton.LeftButton  # set to None to disable
#     viewer.imageViewer.zoomOutButton = Qt.MouseButton.RightButton  # set to None to disable
#     viewer.imageViewer.wheelZoomFactor = 1.25  # Set to None or 1 to disable
#     viewer.imageViewer.panButton = Qt.MouseButton.MiddleButton  # set to None to disable
        
#     # Load an image stack file to be displayed (will popup a file dialog).
#     viewer.open()
        
#     # Show the viewer and run the application.
#     viewer.show()
#     sys.exit(app.exec())
