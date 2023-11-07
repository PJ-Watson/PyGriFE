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
from PyQt6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSlider, QLabel, QComboBox
from PyQt6.QtGui import QPalette, QColor, QImage, QPixmap, QPainter
from QtImageViewer import QtImageViewer
from pathlib import Path
import astropy.io.fits as pf
from astropy.visualization import (LogStretch, MinMaxInterval, PercentileInterval, SqrtStretch, LinearStretch)
import qimage2ndarray
import numpy as np


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

class Color(QWidget):

    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(color))
        self.setPalette(palette)

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("Object Selection")

        layout_h = QHBoxLayout()

        layout_side = QVBoxLayout()
        layout_side.addWidget(Color('green'))
        layout_side.addWidget(Color('blue'))

        # self.direct_interval = MinMaxInterval()

        # Image stretch
        self.stretch = LogStretch()
        stretch_box = QComboBox()
        stretch_box.addItems(["Logarithmic", "Square Root", "Linear"])
        stretch_box.currentTextChanged.connect( self.stretch_update )

        stretch_label = QLabel("Stretch:", self)
        layout_side.addWidget(stretch_label)
        layout_side.addWidget(stretch_box)

        # Interval
        self.interval = PercentileInterval(99.9)
        interval_box = QComboBox()
        interval_box.addItems(["minmax", "99.9%", "99.5%", "99%", "98%", "95%"])
        interval_box.setCurrentText("99.9%")
        interval_box.currentTextChanged.connect( self.interval_update )

        interval_label = QLabel("Interval:", self)
        layout_side.addWidget(interval_label)
        layout_side.addWidget(interval_box)

        # Opacity
        self.opacity = 100
        opacity_box = QComboBox()
        opacity_box.addItems(["100%", "90%", "75%", "50%", "25%", "0%"])
        opacity_box.currentTextChanged.connect( self.opacity_update )

        opacity_label = QLabel("Opacity:", self)
        layout_side.addWidget(opacity_label)
        layout_side.addWidget(opacity_box)

        self.seg_text = QLabel(f"x={0: <6}y={0: <6}\nID={None}", self)
        layout_side.addWidget(self.seg_text)

        self.selected_ids = []
        # # layout_side.addWidget(interval_slider)
        # # layout_side.addWidget(self.interval_label)

        # self.alpha_slider = QSlider(Qt.Orientation.Horizontal, self)
        # self.alpha_slider.setRange(0, 100)
        # self.alpha_slider.setValue(50)
        # self.alpha_slider.setSingleStep(5)
        # self.alpha_slider.setPageStep(10)
        # self.alpha_slider.setTickPosition(QSlider.TickPosition.TicksAbove)

        # self.alpha_slider.valueChanged.connect(self.alpha_update)

        # self.result_label = QLabel('', self)

        # layout_side.addWidget(self.alpha_slider)
        # layout_side.addWidget(self.result_label)
        # # layout_side.setMinimumWidth(200)

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
            if seg_id == 0:
                pass
            if self.selected_ids == [0]:
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

        # seg_id = [1863, 1864]
        # seg_id = 3729

        print (seg_id)
        seg_map = np.isin(self.seg_data, seg_id)#.astype('uint8')
        # seg_map = np.in1d(self.seg_data.ravel(), seg_id).reshape(self.seg_data.shape)#.astype('uint8')
        # seg_map = np.in1d(self.seg_data.ravel(), seg_id, kind="table").reshape(self.seg_data.shape)#.astype('uint8')
        # seg_map = seg_id in self.seg_data
        # print (seg_map)
        # seg_map[self.seg_data!=seg_id] = 0
        # seg_map[self.seg_data==seg_id] = 1
        print (time.time() -t1)

        # seg_plot = self.seg_overlay.copy()
        seg_plot = self.seg_q.copy()
        # seg_plot[:,:,-1] = seg_map
        # seg_plot.setAlphaChannel(qimage2ndarray.array2qimage(seg_map,True))
        seg_plot.setAlphaChannel(qimage2ndarray.gray2qimage(seg_map,True))
        # seg_plot.setAlphaChannel(seg_map)
        # seg_alpha = qimage2ndarray.alpha_view(seg_plot)
        # seg_alpha = seg_map
        # print (dir(seg_alpha))
        # seg_alpha.itemset(seg_alpha)
        # seg_plot2 = qimage2ndarray.array2qimage(seg_plot, True)
        # pixmap = QPixmap.fromImage(seg_plot2)
        pixmap = QPixmap.fromImage(seg_plot)
        # pi

        if not hasattr(self.viewer, "overlay"):
            print ("Creating overlay")
            # qimage = qimage2ndarray.array2qimage(
            #     np.stack(
            #         [
            #             np.zeros_like(seg_map), 
            #             np.ones_like(seg_map), 
            #             np.zeros_like(seg_map), 
            #             seg_map*0.5,
            #         ],
            #         axis=-1), True)
            # pixmap = QPixmap.fromImage(qimage)
            self.viewer.overlay = self.viewer.scene.addPixmap(pixmap)
        else:
            print ("modifying")
            # qimage = qimage2ndarray.array2qimage(
            #     np.stack(
            #         [
            #             np.zeros_like(seg_map), 
            #             np.ones_like(seg_map), 
            #             np.zeros_like(seg_map), 
            #             seg_map*0.5,
            #         ],
            #         axis=-1), True)
            # pixmap = QPixmap.fromImage(qimage)
            self.viewer.overlay.setPixmap(pixmap)

        print (time.time() -t1)

            


    def seg_text_update(self, pos):

        x = int(pos.x())
        y = int(pos.y())
        seg_id = self.seg_data[y, x]

        self.seg_text.setText(f"x={x: <6}y={y: <6}\nID={seg_id}")

        # # print (self.viewer.image())

        # if self.current_seg_id != seg_id:
        #     print (f"Time to update! {self.current_seg_id} -> {seg_id}")
        #     self.current_seg_id = seg_id

        #     painter = QPainter()
        #     painter.begin(self.viewer.image())

        #     # # self.viewer.image().setAlphaChannel(
        #     # # qimage2ndarray.alpha_view(self.viewer.image())*0.9
        #     # # self.alpha_channel[self.seg_data!=seg_id] = 0
        #     seg_map = self.seg_data.copy().astype('uint8')
        #     seg_map[self.seg_data!=seg_id] = 0
        #     seg_map[self.seg_data==seg_id] = 1

        #     try:
        #         painter.drawImage(self.viewer.zoomStack[-1], qimage2ndarray.array2qimage(seg_map))
        #     except Exception as e:
        #         print (e)
        #         painter.drawImage(self.viewer.sceneRect(), qimage2ndarray.array2qimage(seg_map))

        #     painter.end()
            # for x in (self.seg_data==seg_id).nonzero()[0]:
            #     for y in (self.seg_data==seg_id).nonzero()[1]:
            #         print (x,y)
            # r, g, b, a = QColor(self.viewer.image().pixel(x ,y)).getRgb()
            # # ... do something to r, g, b, a ...
            # # print
            # self.viewer.image().setPixel(x, y, QColor(255,255, 255, 0).rgb())
            # self.viewer.repaint()
            # self.update()
            # self.alpha_channel *= 0
            # self.alpha_channel += seg_map*255
            # print (self.alpha_channel)
            # print (qimage2ndarray.alpha_view(self.viewer.image() ))
            # # self.viewer.setImage(self.viewer.image())
            # self.viewer.repaint()
            # self.viewer.update()
            # self.repaint()
            # self.update()
            # # qimage = qimage2ndarray.array2qimage(np.stack([np.full_like(seg_map, 1), seg_map], axis=-1), True)
            # # pixmap = QPixmap.fromImage(qimage)
            # # seg_copy = self.q_seg_img.copy()
            # # seg_copy.setAlphaChannel()
            # # seg_view = qimage2ndarray.recarray_view(seg_copy)
            # # seg_view['a'][self.seg_data!=3729]
            
            # # pixmap = QPixmap.fromImage(seg_copy)

            # if hasattr(self.viewer, "overlay"):
            #     print ("y")
            #     self.viewer.overlay.setPixmap(pixmap)
            # else:
            #     print ("yes")
            #     self.viewer.overlay = self.viewer.scene.addPixmap(pixmap)
            # self.viewer.image().setAlphaChannel(qimage2ndarray.array2qimage(self.alpha_channel))
            # self.repaint()
            # self.update()
            # self.viewer.setImage(self.viewer.image())

        # print (qimage2ndarray.recarray_view(self.viewer.image))


    def interval_update(self, value):
        match value:
            case "minmax":
                self.interval = PercentileInterval(100)
            case _:
                self.interval = PercentileInterval(float(value.split("%")[0]))
        self.reload_image()

    def stretch_update(self, value):
        match value:
            case "Linear":
                self.stretch = LinearStretch()
            case "Square Root":
                self.stretch = SqrtStretch()
            case "Logarithmic":
                self.stretch = LogStretch()
        self.reload_image()

    def opacity_update(self, value):
        self.opacity = float(value.split("%")[0])/100
        self.reload_image()


    def reload_image(self):

        img_array = np.stack(
            [
                self.stretch(self.interval(self.data_array[:,:,0])), 
                np.clip(self.opacity*(1-self.seg_mask)+self.seg_mask, a_min=0, a_max=1)
                #self.data_array[:,:,-1]
            ],
            axis=-1,
        )
        self.viewer.setImage(img_array)
        self.alpha_channel = qimage2ndarray.alpha_view(self.viewer.image())
        print (self.alpha_channel)

    def load_image(self):
        seg_dir = Path("/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3/Prep")
        with pf.open(seg_dir / "nis-wfss-ir_seg.fits") as hdul_seg:
            self.seg_mask = (hdul_seg[0].data > 0).astype(float)
            self.seg_data = hdul_seg[0].data#.astype(int)
            print (self.seg_data.dtype)
            # self.seg_data = self.seg_data [3000:4000,3000:4000]
            # self.seg_mask = self.seg_mask[3000:4000,3000:4000]

            self.seg_overlay = np.stack(
                    [
                        np.zeros_like(self.seg_mask), 
                        np.ones_like(self.seg_mask), 
                        np.zeros_like(self.seg_mask), 
                        self.seg_mask*0.5,
                    ],
                    axis=-1)
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
        with pf.open(seg_dir / "nis-wfss-ir_drz_sci.fits") as hdul_drz:
            print ("here")
            direct_data = hdul_drz[0].data
            # direct_data = direct_data [3000:4000,3000:4000]
            # print (direct_data)
        self.data_array = np.stack([direct_data,self.seg_mask], axis=-1)
        img_array = self.data_array.copy()
        img_array[:,:,0] = self.stretch(self.interval(self.data_array[:,:,0]))
        # print (data_array.shape)

        # self.q_seg_img = qimage2ndarray.array2qimage(np.stack([np.full_like(self.seg_mask, 10000), self.seg_mask], axis=-1), normalize=True,)
        # self.q_seg_view = qimage2ndarray.recarray_view(self.q_seg_img)

        self.q_img = qimage2ndarray.array2qimage(img_array, normalize=True)
        self.viewer.setImage(self.q_img)
        self.q_img_arr = qimage2ndarray.recarray_view(self.q_img)
        self.alpha_channel = qimage2ndarray.alpha_view(self.viewer.image())
        print (self.q_img_arr)



if __name__=="__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()


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
