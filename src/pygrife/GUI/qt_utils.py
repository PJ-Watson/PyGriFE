"""
QtImageViewer.py: PyQt image viewer widget based on QGraphicsView with mouse zooming/panning and ROIs.
Original version by Marcel Goldschen-Ohm, modified by PJW (scroll direction, zoom location, etc).
Also includes worker and signals to allow multithreading.
"""

import os.path
import sys
import traceback
from queue import Queue

from PyQt6.QtCore import (
    QEvent,
    QObject,
    QPoint,
    QPointF,
    QRectF,
    QRunnable,
    QSize,
    Qt,
    QThreadPool,
    pyqtSignal,
    pyqtSlot,
)
from PyQt6.QtGui import (
    QImage,
    QMouseEvent,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QTextCursor,
)
from PyQt6.QtWidgets import (
    QFileDialog,
    QFormLayout,
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsPolygonItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsView,
    QPushButton,
    QSizePolicy,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

# numpy is optional: only needed if you want to display numpy 2d arrays as images.
try:
    import numpy as np
except ImportError:
    np = None

# qimage2ndarray is optional: useful for displaying numpy 2d arrays as images.
# !!! qimage2ndarray requires PyQt5.
#     Some custom code in the viewer appears to handle the conversion from numpy 2d arrays,
#     so qimage2ndarray probably is not needed anymore. I've left it here just in case.
try:
    import qimage2ndarray
except ImportError:
    qimage2ndarray = None

original_author = "Marcel Goldschen-Ohm <marcel.goldschen@gmail.com>"
original_version = "2.0.0"
# Modified in part by PJW


class QtImageViewer(QGraphicsView):
    """
    A PyQt image viewer widget based on QGraphicsView.

    This widget includes mouse zooming/panning and ROIs.

    Class Attributes
    ----------------
    leftMouseButtonPressed : `PyQt6.QtCore.pyqtSignal`
        Test
    """

    # PyQt image viewer widget based on QGraphicsView with mouse zooming/panning and ROIs.

    # Image File:
    # -----------
    # Use the open("path/to/file") method to load an image file into the viewer.
    # Calling open() without a file argument will popup a file selection dialog.

    # Image:
    # ------
    # Use the setImage(im) method to set the image data in the viewer.
    #     - im can be a QImage, QPixmap, or NumPy 2D array (the later requires the package qimage2ndarray).
    #     For display in the QGraphicsView the image will be converted to a QPixmap.

    # Some useful image format conversion utilities:
    #     qimage2ndarray: NumPy ndarray <==> QImage    (https://github.com/hmeine/qimage2ndarray)
    #     ImageQt: PIL Image <==> QImage  (https://github.com/python-pillow/Pillow/blob/master/PIL/ImageQt.py)

    # Mouse:
    # ------
    # Mouse interactions for zooming and panning is fully customizable by simply setting the desired button interactions:
    # e.g.,
    #     regionZoomButton = Qt.LeftButton  # Drag a zoom box.
    #     zoomOutButton = Qt.RightButton  # Pop end of zoom stack (double click clears zoom stack).
    #     panButton = Qt.MiddleButton  # Drag to pan.
    #     wheelZoomFactor = 1.25  # Set to None or 1 to disable mouse wheel zoom.

    # To disable any interaction, just disable its button.
    # e.g., to disable panning:
    #     panButton = None

    # ROIs:
    # -----
    # Can also add ellipse, rectangle, line, and polygon ROIs to the image.
    # ROIs should be derived from the provided EllipseROI, RectROI, LineROI, and PolygonROI classes.
    # ROIs are selectable and optionally moveable with the mouse (see setROIsAreMovable).

    # TODO: Add support for editing the displayed image contrast.
    # TODO: Add support for drawing ROIs with the mouse.
    #

    # Mouse button signals emit image scene (x, y) coordinates.
    # > For image (row, column) matrix indexing, row = y and column = x.
    # > These signals will NOT be emitted if the event is handled by an interaction such as zoom or pan.
    # > If aspect ratio prevents image from filling viewport, emitted position may be outside image bounds.
    leftMouseButtonPressed = pyqtSignal(float, float)
    leftMouseButtonReleased = pyqtSignal(Qt.KeyboardModifier, float, float)
    middleMouseButtonPressed = pyqtSignal(float, float)
    middleMouseButtonReleased = pyqtSignal(float, float)
    rightMouseButtonPressed = pyqtSignal(float, float)
    rightMouseButtonReleased = pyqtSignal(float, float)
    leftMouseButtonDoubleClicked = pyqtSignal(float, float)
    rightMouseButtonDoubleClicked = pyqtSignal(float, float)

    # Emitted upon zooming/panning.
    viewChanged = pyqtSignal()

    # Emitted on mouse motion.
    # Emits mouse position over image in image pixel coordinates.
    # > setMouseTracking(True) if you want to use this at all times.
    mousePositionOnImageChanged = pyqtSignal(QPoint)

    # Emit index of selected ROI
    roiSelected = pyqtSignal(int)

    def __init__(self):
        QGraphicsView.__init__(self)

        # Image is displayed as a QPixmap in a QGraphicsScene attached to this QGraphicsView.
        self.scene = QGraphicsScene()
        self.setScene(self.scene)

        # Better quality pixmap scaling?
        # self.setRenderHints(QPainter.Antialiasing | QPainter.SmoothPixmapTransform)

        # Displayed image pixmap in the QGraphicsScene.
        self._image = None

        # Image aspect ratio mode.
        #   Qt.IgnoreAspectRatio: Scale image to fit viewport.
        #   Qt.KeepAspectRatio: Scale image to fit inside viewport, preserving aspect ratio.
        #   Qt.KeepAspectRatioByExpanding: Scale image to fill the viewport, preserving aspect ratio.
        self.aspectRatioMode = Qt.AspectRatioMode.KeepAspectRatio

        # Scroll bar behaviour.
        #   Qt.ScrollBarAlwaysOff: Never shows a scroll bar.
        #   Qt.ScrollBarAlwaysOn: Always shows a scroll bar.
        #   Qt.ScrollBarAsNeeded: Shows a scroll bar only when zoomed.
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        # Interactions (set buttons to None to disable interactions)
        # > Events handled by interactions will NOT emit `MouseButton`` signals.
        # > Note: regionZoomButton will still emit a `MouseButtonReleased`` signal on a click (i.e. tiny box).
        self.regionZoomButton = Qt.MouseButton.LeftButton  # Drag a zoom box.
        self.zoomOutButton = (
            Qt.MouseButton.RightButton
        )  # Pop end of zoom stack (double click clears zoom stack).
        self.panButton = Qt.MouseButton.MiddleButton  # Drag to pan.
        self.wheelZoomFactor = 1.25  # Set to None or 1 to disable mouse wheel zoom.

        # Stack of QRectF zoom boxes in scene coordinates.
        # > If you update this manually, be sure to call updateViewer() to reflect any changes.
        self.zoomStack = []

        # Flags for active zooming/panning.
        self._isZooming = False
        self._isPanning = False

        # Store temporary position in screen pixels or scene units.
        self._pixelPosition = QPoint()
        self._scenePosition = QPointF()

        # Track mouse position. e.g., For displaying coordinates in a UI.
        # self.setMouseTracking(True)

        # ROIs.
        self.ROIs = []

        # # For drawing ROIs.
        # self.drawROI = None

        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def sizeHint(self):
        """
        Returns a hint for the size of the current scene.

        Returns
        -------
        `PyQt6.QtCore.QSize`
            A suggested size.
        """
        return QSize(900, 600)

    def hasImage(self):
        """
        Check if the scene contains an image pixmap.

        Returns
        -------
        bool
            ``True`` if the scene contains an image.
        """
        return self._image is not None

    def clearImage(self):
        """
        Removes the current image pixmap from the scene if it exists.
        """
        if self.hasImage():
            self.scene.removeItem(self._image)
            self._image = None

    def pixmap(self):
        """
        Returns the scene's current image pixmap.

        Returns
        -------
        `PyQt6.QtGui.QPixmap` or None
            The current scene image, if it exists.
        """
        if self.hasImage():
            return self._image.pixmap()
        return None

    def image(self):
        """
        Returns the scene's current image as a QImage.

        Returns
        -------
        `PyQt6.QtGui.QImage` or None
            The current scene image, if it exists.
        """
        if self.hasImage():
            return self._image.pixmap().toImage()
        return None

    def setImage(self, image):
        """
        Set the scene's current image pixmap to the input.

        Parameters
        ----------
        image :  `PyQt.QtGui.QImage` or `PyQt.QtGui.QPixmap`
            The input image.

        Raises
        ------
        RuntimeError
            Raised if the input image has type other than
            `PyQt6.QtGui.QImage` or `PyQt6.QtGui.QPixmap`.
        """

        if type(image) is QPixmap:
            pixmap = image
        elif type(image) is QImage:
            pixmap = QPixmap.fromImage(image)
        elif (np is not None) and (type(image) is np.ndarray):
            if qimage2ndarray is not None:
                qimage = qimage2ndarray.array2qimage(image, True)
                pixmap = QPixmap.fromImage(qimage)
            else:
                image = image.astype(np.float32)
                image -= image.min()
                image /= image.max()
                image *= 255
                image[image > 255] = 255
                image[image < 0] = 0
                image = image.astype(np.uint8)
                height, width = image.shape
                bytes = image.tobytes()
                qimage = QImage(bytes, width, height, QImage.Format.Format_Grayscale8)
                pixmap = QPixmap.fromImage(qimage)
        else:
            raise RuntimeError(
                "ImageViewer.setImage: Argument must be a QImage, QPixmap, or"
                " numpy.ndarray."
            )
        if self.hasImage():
            self._image.setPixmap(pixmap)
        else:
            self._image = self.scene.addPixmap(pixmap)

        # Better quality pixmap scaling?
        # !!! This will distort actual pixel data when zoomed way in.
        #     For scientific image analysis, you probably don't want this.
        # self._pixmap.setTransformationMode(Qt.SmoothTransformation)

        self.setSceneRect(QRectF(pixmap.rect()))  # Set scene size to image size.
        self.updateViewer()

    def open(self, filepath=None):
        """
        Load an image from file.

        Parameters
        ----------
        filepath : str or os.PathLike, optional
            The filepath of the image to load. If no image is supplied,
            a file dialogue will be opened to choose the image file.
        """
        if filepath is None:
            filepath, dummy = QFileDialog.getOpenFileName(self, "Open image file.")
        if len(filepath) and os.path.isfile(filepath):
            image = QImage(filepath)
            self.setImage(image)

    def updateViewer(self):
        """
        Show the current zoom.

        If the entire image is visible, apply current aspect ratio mode.
        """
        if not self.hasImage():
            return
        if len(self.zoomStack):
            self.fitInView(
                self.zoomStack[-1], self.aspectRatioMode
            )  # Show zoomed rect.
        else:
            self.fitInView(self.sceneRect(), self.aspectRatioMode)  # Show entire image.

    def clearZoom(self):
        """
        Clear the zoom status of the scene.
        """
        if len(self.zoomStack) > 0:
            self.zoomStack = []
            self.updateViewer()
            self.viewChanged.emit()

    def resizeEvent(self, event):
        """
        Maintain the current zoom on resize.
        """
        self.updateViewer()

    def mousePressEvent(self, event):
        """
        Start either the mouse pan or zoom mode.

        Parameters
        ----------
        event : `PyQt6.QtGui.QMouseEvent`
            The event that triggered this, i.e. a mouse press.
        """

        # Ignore dummy events. e.g., Faking pan with left button ScrollHandDrag.
        dummyModifiers = Qt.KeyboardModifier(
            Qt.KeyboardModifier.ShiftModifier
            | Qt.KeyboardModifier.ControlModifier
            | Qt.KeyboardModifier.AltModifier
            | Qt.KeyboardModifier.MetaModifier
        )
        if event.modifiers() == dummyModifiers:
            QGraphicsView.mousePressEvent(self, event)
            event.accept()
            return

        # # Draw ROI
        # if self.drawROI is not None:
        #     if self.drawROI == "Ellipse":
        #         # Click and drag to draw ellipse. +Shift for circle.
        #         pass
        #     elif self.drawROI == "Rect":
        #         # Click and drag to draw rectangle. +Shift for square.
        #         pass
        #     elif self.drawROI == "Line":
        #         # Click and drag to draw line.
        #         pass
        #     elif self.drawROI == "Polygon":
        #         # Click to add points to polygon. Double-click to close polygon.
        #         pass

        # Start dragging a region zoom box?
        if (self.regionZoomButton is not None) and (
            event.button() == self.regionZoomButton
        ):
            self._pixelPosition = event.pos()  # store pixel position
            self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
            QGraphicsView.mousePressEvent(self, event)
            event.accept()
            self._isZooming = True
            return

        if (self.zoomOutButton is not None) and (event.button() == self.zoomOutButton):
            if len(self.zoomStack):
                self.zoomStack.pop()
                self.updateViewer()
                self.viewChanged.emit()
            event.accept()
            return

        # Start dragging to pan?
        if (self.panButton is not None) and (event.button() == self.panButton):
            self._pixelPosition = event.pos()  # store pixel position
            self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
            if self.panButton == Qt.MouseButton.LeftButton:
                QGraphicsView.mousePressEvent(self, event)
            else:
                # ScrollHandDrag ONLY works with LeftButton, so fake it.
                # Use a bunch of dummy modifiers to notify that event should NOT be handled as usual.
                self.viewport().setCursor(Qt.CursorShape.ClosedHandCursor)
                dummyModifiers = Qt.KeyboardModifier(
                    Qt.KeyboardModifier.ShiftModifier
                    | Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.AltModifier
                    | Qt.KeyboardModifier.MetaModifier
                )
                dummyEvent = QMouseEvent(
                    QEvent.Type.MouseButtonPress,
                    QPointF(event.pos()),
                    Qt.MouseButton.LeftButton,
                    event.buttons(),
                    dummyModifiers,
                )
                self.mousePressEvent(dummyEvent)
            sceneViewport = (
                self.mapToScene(self.viewport().rect())
                .boundingRect()
                .intersected(self.sceneRect())
            )
            self._scenePosition = sceneViewport.topLeft()
            event.accept()
            self._isPanning = True
            return

        scenePos = self.mapToScene(event.pos())
        if event.button() == Qt.MouseButton.LeftButton:
            self.leftMouseButtonPressed.emit(scenePos.x(), scenePos.y())
        elif event.button() == Qt.MouseButton.MiddleButton:
            self.middleMouseButtonPressed.emit(scenePos.x(), scenePos.y())
        elif event.button() == Qt.MouseButton.RightButton:
            self.rightMouseButtonPressed.emit(scenePos.x(), scenePos.y())

        QGraphicsView.mousePressEvent(self, event)

    def mouseReleaseEvent(self, event):
        """
        Stop mouse pan or zoom mode, and apply zoom if valid.

        Parameters
        ----------
        event : `PyQt6.QtGui.QMouseEvent`
            The event that triggered this, i.e. a mouse release.
        """

        # Ignore dummy events. e.g., Faking pan with left button ScrollHandDrag.
        dummyModifiers = Qt.KeyboardModifier(
            Qt.KeyboardModifier.ShiftModifier
            | Qt.KeyboardModifier.ControlModifier
            | Qt.KeyboardModifier.AltModifier
            | Qt.KeyboardModifier.MetaModifier
        )
        if event.modifiers() == dummyModifiers:
            QGraphicsView.mouseReleaseEvent(self, event)
            event.accept()
            return

        # Finish dragging a region zoom box?
        if (self.regionZoomButton is not None) and (
            event.button() == self.regionZoomButton
        ):
            QGraphicsView.mouseReleaseEvent(self, event)
            zoomRect = (
                self.scene.selectionArea().boundingRect().intersected(self.sceneRect())
            )
            # Clear current selection area (i.e. rubberband rect).
            self.scene.setSelectionArea(QPainterPath())
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
            # If zoom box is 3x3 screen pixels or smaller, do not zoom and proceed to process as a click release.
            zoomPixelWidth = abs(event.pos().x() - self._pixelPosition.x())
            zoomPixelHeight = abs(event.pos().y() - self._pixelPosition.y())
            if zoomPixelWidth > 3 and zoomPixelHeight > 3:
                if zoomRect.isValid() and (zoomRect != self.sceneRect()):
                    self.zoomStack.append(zoomRect)
                    self.updateViewer()
                    self.viewChanged.emit()
                    event.accept()
                    self._isZooming = False
                    return

        # Finish panning?
        if (self.panButton is not None) and (event.button() == self.panButton):
            if self.panButton == Qt.MouseButton.LeftButton:
                QGraphicsView.mouseReleaseEvent(self, event)
            else:
                # ScrollHandDrag ONLY works with LeftButton, so fake it.
                # Use a bunch of dummy modifiers to notify that event should NOT be handled as usual.
                self.viewport().setCursor(Qt.CursorShape.ArrowCursor)
                dummyModifiers = Qt.KeyboardModifier(
                    Qt.KeyboardModifier.ShiftModifier
                    | Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.AltModifier
                    | Qt.KeyboardModifier.MetaModifier
                )
                dummyEvent = QMouseEvent(
                    QEvent.Type.MouseButtonRelease,
                    QPointF(event.pos()),
                    Qt.MouseButton.LeftButton,
                    event.buttons(),
                    dummyModifiers,
                )
                self.mouseReleaseEvent(dummyEvent)
            self.setDragMode(QGraphicsView.DragMode.NoDrag)
            if len(self.zoomStack) > 0:
                sceneViewport = (
                    self.mapToScene(self.viewport().rect())
                    .boundingRect()
                    .intersected(self.sceneRect())
                )
                delta = sceneViewport.topLeft() - self._scenePosition
                self.zoomStack[-1].translate(delta)
                self.zoomStack[-1] = self.zoomStack[-1].intersected(self.sceneRect())
                self.viewChanged.emit()
            event.accept()
            self._isPanning = False
            return

        scenePos = self.mapToScene(event.pos())
        if event.button() == Qt.MouseButton.LeftButton:
            self.leftMouseButtonReleased.emit(
                event.modifiers(), scenePos.x(), scenePos.y()
            )
        elif event.button() == Qt.MouseButton.MiddleButton:
            self.middleMouseButtonReleased.emit(scenePos.x(), scenePos.y())
        elif event.button() == Qt.MouseButton.RightButton:
            self.rightMouseButtonReleased.emit(scenePos.x(), scenePos.y())

        QGraphicsView.mouseReleaseEvent(self, event)

    def mouseDoubleClickEvent(self, event):
        """
        Reset the zoom to show the entire image.

        Parameters
        ----------
        event : `PyQt6.QtGui.QMouseEvent`
            The event that triggered this, i.e. a double click.
        """
        # Zoom out on double click?
        if (self.zoomOutButton is not None) and (event.button() == self.zoomOutButton):
            self.clearZoom()
            event.accept()
            return

        scenePos = self.mapToScene(event.pos())
        if event.button() == Qt.MouseButton.LeftButton:
            self.leftMouseButtonDoubleClicked.emit(scenePos.x(), scenePos.y())
        elif event.button() == Qt.MouseButton.RightButton:
            self.rightMouseButtonDoubleClicked.emit(scenePos.x(), scenePos.y())

        QGraphicsView.mouseDoubleClickEvent(self, event)

    def wheelEvent(self, event):
        """
        Handle all mouse wheel scroll events.

        If the shift key is pressed, zoom into/out of the image; else
        scroll around the image in the current zoom state.

        Parameters
        ----------
        event : `PyQt6.QtGui.QWheelEvent`
            The event that triggered this, i.e. a mouse wheel scroll.
        """
        if (
            self.wheelZoomFactor is not None
            and event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        ):
            if self.wheelZoomFactor == 1:
                return
            if event.angleDelta().y() > 0:
                # zoom in
                if len(self.zoomStack) == 0:
                    self.zoomStack.append(self.sceneRect())
                elif len(self.zoomStack) > 1:
                    del self.zoomStack[:-1]
                zoomRect = self.zoomStack[-1]
                center = zoomRect.center()
                new_center = self.mapToScene(event.position().toPoint())  # Scene coords
                zoomRect.setLeft(
                    new_center.x()
                    - (new_center.x() - zoomRect.left()) / self.wheelZoomFactor
                )
                zoomRect.setRight(
                    new_center.x()
                    + (zoomRect.right() - new_center.x()) / self.wheelZoomFactor
                )
                zoomRect.setBottom(
                    new_center.y()
                    - (new_center.y() - zoomRect.bottom()) / self.wheelZoomFactor
                )
                zoomRect.setTop(
                    new_center.y()
                    + (zoomRect.top() - new_center.y()) / self.wheelZoomFactor
                )

                self.zoomStack[-1] = zoomRect.intersected(self.sceneRect())
                self.updateViewer()
                self.viewChanged.emit()
            else:
                # zoom out
                if len(self.zoomStack) == 0:
                    # Already fully zoomed out.
                    return
                if len(self.zoomStack) > 1:
                    del self.zoomStack[:-1]
                zoomRect = self.zoomStack[-1]
                center = zoomRect.center()

                new_center = self.mapToScene(event.position().toPoint())  # Scene coords
                zoomRect.setLeft(
                    new_center.x()
                    - (new_center.x() - zoomRect.left()) * self.wheelZoomFactor
                )
                zoomRect.setRight(
                    new_center.x()
                    + (zoomRect.right() - new_center.x()) * self.wheelZoomFactor
                )
                zoomRect.setBottom(
                    new_center.y()
                    - (new_center.y() - zoomRect.bottom()) * self.wheelZoomFactor
                )
                zoomRect.setTop(
                    new_center.y()
                    + (zoomRect.top() - new_center.y()) * self.wheelZoomFactor
                )

                self.zoomStack[-1] = zoomRect.intersected(self.sceneRect())
                if self.zoomStack[-1] == self.sceneRect():
                    self.zoomStack = []
                self.updateViewer()
                self.viewChanged.emit()
            event.accept()
            return

        QGraphicsView.wheelEvent(self, event)

    def mouseMoveEvent(self, event):
        """
        Handle all mouse movement in the scene.

        Emit the mouse position, and change the view if panning.

        Parameters
        ----------
        event : `PyQt6.QtGui.QMouseEvent`
            The mouse move event.
        """
        # Emit updated view during panning.
        if self._isPanning:
            QGraphicsView.mouseMoveEvent(self, event)
            if len(self.zoomStack) > 0:
                sceneViewport = (
                    self.mapToScene(self.viewport().rect())
                    .boundingRect()
                    .intersected(self.sceneRect())
                )
                delta = sceneViewport.topLeft() - self._scenePosition
                self._scenePosition = sceneViewport.topLeft()
                self.zoomStack[-1].translate(delta)
                self.zoomStack[-1] = self.zoomStack[-1].intersected(self.sceneRect())
                self.updateViewer()
                self.viewChanged.emit()

        scenePos = self.mapToScene(event.pos())
        if self.sceneRect().contains(scenePos):
            # Pixel index offset from pixel center.
            x = int(round(scenePos.x() - 0.5))
            y = int(round(scenePos.y() - 0.5))
            imagePos = QPoint(x, y)
        else:
            # Invalid pixel position.
            imagePos = QPoint(-1, -1)
        self.mousePositionOnImageChanged.emit(imagePos)

        QGraphicsView.mouseMoveEvent(self, event)

    def enterEvent(self, event):
        """
        Change the cursor shape when the mouse enters the scene.

        Parameters
        ----------
        event : `PyQt6.QtGui.QEnterEvent`
            The triggering event.
        """
        self.setCursor(Qt.CursorShape.CrossCursor)

    def leaveEvent(self, event):
        """
        Change the cursor shape when the mouse leaves the scene.

        Parameters
        ----------
        event : `PyQt6.QtCore.QEvent`
            The triggering event.
        """
        self.setCursor(Qt.CursorShape.ArrowCursor)

    def addROIs(self, rois):
        """
        Add ROIs to the scene.

        Parameters
        ----------
        rois : `PyQt6.QtWidgets.QtGraphicsItem`
            The regions of interest to add to the scene.
        """
        for roi in rois:
            self.scene.addItem(roi)
            self.ROIs.append(roi)

    def deleteROIs(self, rois):
        """
        Remove ROIs from the scene.

        Parameters
        ----------
        rois : `PyQt6.QtWidgets.QtGraphicsItem`
            The regions of interest to remove from the scene.
        """
        for roi in rois:
            self.scene.removeItem(roi)
            self.ROIs.remove(roi)
            del roi

    def clearROIs(self):
        """
        Remove all ROIs from the scene.
        """
        for roi in self.ROIs:
            self.scene.removeItem(roi)
        del self.ROIs[:]

    def roiClicked(self, roi):
        """
        Emit a signal when a ROI is clicked.

        Parameters
        ----------
        roi : `PyQt6.QtWidgets.QtGraphicsItem`
            The selected ROI.
        """
        for i in range(len(self.ROIs)):
            if roi is self.ROIs[i]:
                self.roiSelected.emit(i)
                print(i)
                break

    def setROIsAreMovable(self, tf):
        """
        Set movable ROIs.

        Parameters
        ----------
        tf : bool
            If the ROIs are movable.
        """
        if tf:
            for roi in self.ROIs:
                roi.setFlags(roi.flags() | QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        else:
            for roi in self.ROIs:
                roi.setFlags(
                    roi.flags() & ~QGraphicsItem.GraphicsItemFlag.ItemIsMovable
                )

    def addSpots(self, xy, radius):
        """
        Add circular ROIs at the given coordinates.

        Parameters
        ----------
        xy : array-like
            A set of coordinates.
        radius : float
            The radius of the ROI.
        """
        for xy_ in xy:
            x, y = xy_
            spot = EllipseROI(self)
            spot.setRect(x - radius, y - radius, 2 * radius, 2 * radius)
            self.scene.addItem(spot)
            self.ROIs.append(spot)


class EllipseROI(QGraphicsEllipseItem):
    def __init__(self, viewer):
        QGraphicsItem.__init__(self)
        self._viewer = viewer
        pen = QPen(Qt.yellow)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setFlags(self.GraphicsItemFlag.ItemIsSelectable)

    def mousePressEvent(self, event):
        QGraphicsItem.mousePressEvent(self, event)
        if event.button() == Qt.MouseButton.LeftButton:
            self._viewer.roiClicked(self)


class RectROI(QGraphicsRectItem):
    def __init__(self, viewer):
        QGraphicsItem.__init__(self)
        self._viewer = viewer
        pen = QPen(Qt.GlobalColor.yellow)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setFlags(self.GraphicsItemFlag.ItemIsSelectable)

    def mousePressEvent(self, event):
        QGraphicsItem.mousePressEvent(self, event)
        if event.button() == Qt.MouseButton.LeftButton:
            self._viewer.roiClicked(self)


class LineROI(QGraphicsLineItem):
    def __init__(self, viewer):
        QGraphicsItem.__init__(self)
        self._viewer = viewer
        pen = QPen(Qt.GlobalColor.yellow)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setFlags(self.GraphicsItemFlag.ItemIsSelectable)

    def mousePressEvent(self, event):
        QGraphicsItem.mousePressEvent(self, event)
        if event.button() == Qt.MouseButton.LeftButton:
            self._viewer.roiClicked(self)


class PolygonROI(QGraphicsPolygonItem):
    def __init__(self, viewer):
        QGraphicsItem.__init__(self)
        self._viewer = viewer
        pen = QPen(Qt.GlobalColor.yellow)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setFlags(self.GraphicsItemFlag.ItemIsSelectable)

    def mousePressEvent(self, event):
        QGraphicsItem.mousePressEvent(self, event)
        if event.button() == Qt.MouseButton.LeftButton:
            self._viewer.roiClicked(self)


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    """

    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int, str)


class Worker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    """

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs["progress_callback"] = self.signals.progress

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done


# The new Stream Object which replaces the default stream associated with sys.stdout
# This object just puts data in a queue!
class WriteStream(object):
    def __init__(self, queue):
        self.queue = queue

    def write(self, text):
        self.queue.put(text)


def receiver_fn(queue, progress_callback=None):
    while True:
        text = queue.get()
        progress_callback.emit(text)


class TerminalWindow(QWidget):
    def __init__(self, root):
        super().__init__()
        self.root = root
        self.setWindowTitle("Terminal Output")

        self.v_layout = QVBoxLayout()

        self.textedit = QTextEdit()
        self.v_layout.addWidget(self.textedit)

        self.setLayout(self.v_layout)
        self.setMinimumWidth(720)
        self.setMinimumWidth(568)

    @pyqtSlot(str)
    def append_text(self, text):
        self.textedit.moveCursor(QTextCursor.MoveOperation.End)
        self.textedit.insertPlainText(text)
