
import os
import glob
import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import matplotlib.pyplot as plt
from grizli import utils, prep, jwst_utils, multifit
from grizli.pipeline import auto_script
from grizli.pipeline.auto_script import get_yml_parameters
import grizli
from importlib import reload
from grizli import model, multifit, grismconf
from grizli import fitting
import astropy.units as u
from tqdm import tqdm
from pathlib import Path

"""
Steps:
 - Recalculate photometry for affected objects - look at this function:

        auto_script.multiband_catalog(field_root=root,
            **multiband_catalog_args)

 - Edit contamination maps - need to modify `multifit.GroupFLT`

 - Permanently modify *GrismFLT.fits, or re-run modifications each time?

 - Force all functions here to save to a different folder, e.g. ForcedExtractions, and save all files there.

 - ImageViewer module can be as generic as possible, but `Extract Object` button can be customised to work only with grizli.

 - Sub-class imageviewer? Make one completely generic function, and one specifically for grizli.

"""

def load_contamination_maps(grism_dir, cat):

    grism_dir = Path(grism_dir)

    grp = multifit.GroupFLT(
        grism_files=[*grism_dir.glob("*GrismFLT.fits")],
        catalog="{0}-ir.cat.fits".format(root),
        cpu_count=-1,
        sci_extn=1,
        pad=800,
    )


print("Grizli version: ", grizli.__version__)

print("1. Importing Python packages...[COMPLETE]")
############################################
### Define necessary variables and paths ###
############################################
print("2. Defining variables and paths...")

HOME_PATH = "/media/sharedData/data/GLASS_owncloud/NIRISS/ABELL2744/v3"  # change this
root = "nis-wfss"

os.chdir(HOME_PATH + "/Prep")

print("2. Defining variables and paths...[COMPLETE]")