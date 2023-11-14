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
from seg_map_viewer import SegMapViewer, Separator, FilesWindow, cQLineEdit
import json

class GrizliGUI(SegMapViewer):
    def __init__(self, new_directory="ForcedExtractions", filters=["F115W", "F150W", "F200W"]):
        super().__init__()

        self.layout_side.addWidget(Separator())
        extract_object_button = QPushButton('Extract Object', self)
        extract_object_button.clicked.connect(self.extract_object)
        self.layout_side.addWidget(extract_object_button)

        self.prep_dir = None
        self.field_name = ""
        self.new_directory = new_directory

    def open_files_window(self, event=None):

        if self.files_window is None:
            self.files_window = GrizliFilesWindow(self)
        self.files_window.show()

    def save_output(self, event=None):

        if not (hasattr(self, "seg_img_path") and hasattr(self, "seg_data")):
            print ("No segmentation mask loaded.")
            return

        self.new_dir_path = Path(self.prep_dir.parent) / self.new_directory
        self.new_dir_path.mkdir(exist_ok=True, parents=True)

        with open(self.new_dir_path / "remapped_ids.json", "w") as f:
            json.dump(self.remapped_ids, f)

        with pf.open(self.seg_img_path) as seg_hdul:
            
            seg_hdul[0].data = self.seg_data[::-1,:]

            seg_hdul.writeto(self.new_dir_path / self.seg_img_path.name, overwrite=True)

    def extract_object(self, event=None):

        print (self.selected_ids)


class GrizliFilesWindow(FilesWindow):
    def __init__(self, root):
        super().__init__(root)

        self.prep_dir_line = cQLineEdit(parent=self, is_dir=True)
        self.sub_layout.insertRow(0, "Prep Directory", self.prep_dir_line)

    def change_directory(self, event=None):

        if self.prep_dir_line.text() is None or self.prep_dir_line.text()=="":
            if self.root.prep_dir is not None:
                init = str(prep_dir)
            elif self.recent_dir is None:
                init = str(self.recent_dir)
            else:
                init = str(Path.home())

            dir_name = QFileDialog.getExistingDirectory(self, "Open directory", init)
            if dir_name:
                self.root.prep_dir = Path(dir_name)
                self.prep_dir_line.setText(str(self.root.prep_dir))
            else:
                return
        else:
            self.root.prep_dir = Path(self.prep_dir_line.text())

        try:
            self.seg_line.setText(str([*self.root.prep_dir.glob("*ir_seg.fits")][0]))
            self.root.field_name =  ([*self.root.prep_dir.glob("*ir_seg.fits")][0].stem.split("-ir_seg")[0])
        except:
            print ("Segmentation map not found.")

        try:
            self.stack_line.setText(str([*self.root.prep_dir.glob(f"*{self.root.field_name}-ir_drz_sci.fits")][0]))
        except:
            print ("Stacked image not found.")

        try:
            for f, l in zip(self.root.filters, [self.b_line, self.g_line, self.r_line]):
                l.setText(str([*self.root.prep_dir.glob(f"*{f.lower()}_drz_sci.fits")][0]))
        except:
            print ("Could not find all filter images.")

    def load_all(self):
        self.root.prep_dir = Path(self.prep_dir_line.text())
        super().load_all()
    
if __name__=="__main__":
    app = QApplication(sys.argv)
    window = GrizliGUI()
    window.showMaximized()
    app.exec()