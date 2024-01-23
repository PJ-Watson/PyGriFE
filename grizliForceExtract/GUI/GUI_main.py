import json
import sys
import time

# from .QtImageViewer import QtImageViewer
from pathlib import Path
from queue import Queue

import astropy.io.fits as pf
import numpy as np
import qimage2ndarray
from astropy import wcs
from astropy.visualization import (
    AsinhStretch,
    LinearStretch,
    LogStretch,
    ManualInterval,
    SqrtStretch,
)
from grizli_extractor import GrizliExtractor
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QImage, QPainter, QPalette, QPixmap
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)
from qt_utils import TerminalWindow, Worker, WriteStream
from seg_map_viewer import FilesWindow, SegMapViewer, Separator, cQLineEdit


class GrizliGUI(SegMapViewer):
    def __init__(
        self, new_directory="ForcedExtractions", filters=["F115W", "F150W", "F200W"]
    ):
        super().__init__()

        self.layout_side.addWidget(Separator())
        self.extract_object_button = QPushButton("Extract Object", self)
        self.extract_object_button.clicked.connect(self.extraction_handler)
        self.layout_side.addWidget(self.extract_object_button)

        self.terminal_window = None
        self.extract_in_progress = False

        self.prep_dir = None
        self.field_name = ""
        self.new_directory = new_directory
        self.ge = None

    def open_files_window(self, event=None):
        if self.files_window is None:
            self.files_window = GrizliFilesWindow(self)
        self.files_window.show()

    def open_terminal_window(self, event=None):
        if self.terminal_window is None:
            self.terminal_window = TerminalWindow(self)
        self.terminal_window.show()

    def save_output(self, event=None):
        if not (hasattr(self, "seg_img_path") and hasattr(self, "seg_data")):
            print("No segmentation mask loaded.")
            return

        self.new_dir_path = Path(self.prep_dir.parent) / self.new_directory
        self.new_dir_path.mkdir(exist_ok=True, parents=True)

        with open(self.new_dir_path / "remapped_ids.json", "w") as f:
            json.dump(self.remapped_ids, f)

        with pf.open(self.seg_img_path) as seg_hdul:
            seg_hdul[0].data = self.seg_data[::-1, :]

            seg_hdul.writeto(self.new_dir_path / self.seg_img_path.name, overwrite=True)

    def receiver_fn(self, queue, progress_callback=None):
        while self.extract_in_progress:
            text = queue.get()
            progress_callback.emit(text)
        return

    def extraction_handler(self, event=None):
        self.extract_in_progress = True
        self.extract_object_button.setEnabled(False)
        self.open_terminal_window()

        queue = Queue()
        sys.stdout = WriteStream(queue)

        receive_worker = Worker(self.receiver_fn, queue)
        receive_worker.signals.progress.connect(self.terminal_window.append_text)
        self.threadpool.start(receive_worker)

        extract_worker = Worker(self.extract_object)
        # worker.signals.progress.connect(self.progress_fn)
        # worker.signals.result.connect(self.root.set_img)
        extract_worker.signals.finished.connect(self.finish_extractions)
        self.threadpool.start(extract_worker)

        # self.threadpool

    def finish_extractions(self):
        # print ("pls end")
        self.extract_in_progress = False
        self.extract_object_button.setEnabled(True)
        sys.stdout = sys.__stdout__

    def extract_object(self, event=None, progress_callback=None):
        # import logging
        # root = logging.getLogger()
        # root.setLevel(logging.INFO)
        # fh = logging.FileHandler('debug.log')
        # fh.setLevel(logging.INFO)

        # old_stdout = sys.stdout    # in case you want to restore later
        # sys.stdout = fh.stream

        # root.addHandler(fh)
        print("Beginning extraction.")
        print(self.selected_ids)

        self.new_dir_path = Path(self.prep_dir.parent) / self.new_directory
        self.new_dir_path.mkdir(exist_ok=True, parents=True)

        if self.ge is None:
            self.ge = GrizliExtractor(self.field_name, self.prep_dir, self.new_dir_path)
        self.ge.load_seg_img(self.seg_data[::-1, :])
        self.ge.regen_multiband_catalogue()
        if not hasattr(self.ge, "grp"):
            self.ge.load_contamination_maps()
        # self.ge.extract_sep()
        self.ge.extract_spectra(self.selected_ids)
        # sys.stdout = old_stdout
        return


class GrizliFilesWindow(FilesWindow):
    def __init__(self, root):
        super().__init__(root)

        self.prep_dir_line = cQLineEdit(parent=self, is_dir=True)
        self.sub_layout.insertRow(0, "Prep Directory", self.prep_dir_line)

    def change_directory(self, event=None):
        if self.prep_dir_line.text() is None or self.prep_dir_line.text() == "":
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
            self.root.field_name = [*self.root.prep_dir.glob("*ir_seg.fits")][
                0
            ].stem.split("-ir_seg")[0]
        except:
            print("Segmentation map not found.")

        try:
            self.stack_line.setText(
                str(
                    [
                        *self.root.prep_dir.glob(
                            f"*{self.root.field_name}-ir_drz_sci.fits"
                        )
                    ][0]
                )
            )
        except:
            print("Stacked image not found.")

        try:
            for f, l in zip(self.root.filters, [self.b_line, self.g_line, self.r_line]):
                l.setText(
                    str([*self.root.prep_dir.glob(f"*{f.lower()}_drz_sci.fits")][0])
                )
        except:
            print("Could not find all filter images.")

    def load_all(self):
        self.root.prep_dir = Path(self.prep_dir_line.text())
        super().load_all()
