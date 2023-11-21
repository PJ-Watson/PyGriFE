from setuptools import Extension, setup
import os
import re
import sys
from glob import glob

fname = "grizliForceExtract/utils_c/sep_fns.pyx"
import numpy
from Cython.Build import cythonize
sourcefiles = [fname] + glob(os.path.join("grizliForceExtract/utils_c/sep_fns", "*.c"))
headerfiles = glob(os.path.join("grizliForceExtract/utils_c/sep_fns", "*.h"))
include_dirs = [numpy.get_include(), "grizliForceExtract/utils_c/sep_fns"]
extensions = [
    Extension(
        name="grizliForceExtract.utils_c.sep_fns", 
        sources=sourcefiles,
        include_dirs=include_dirs,
        depends=headerfiles, 
        define_macros=[("_USE_MATH_DEFINES", "1")],
    )
]
extensions = cythonize(extensions)

setup(
    # name="grizliForceExtract",
    ext_modules=extensions,
)