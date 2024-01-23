# PyGriFE

A package to aid in the forced extraction of pre-processed grism data.

## Installation

In all cases, it is strongly recommended to install `PyGriFE` into a new virtual environment, to minimise dependency conflicts (see [Requirements](#requirements)).

### Building from source

To clone the latest GitHub repository, use this command:

```
git clone https://github.com/PJ-Watson/PyGriFE.git
```

To build and install `PyGriFE`, run (from the root of the source tree):

```
pip install .
```

## Necessary Files

`PyGriFE` is designed to work with pre-processed grism data, specifically the output from [`grizli`](https://github.com/gbrammer/grizli). If you do not have the following files, `PyGriFe` will not run:
 - "_*GrismFLT.fits_" and "_*GrismFLT.pkl_" files, containing the grism exposures, contamination models, blotted segmentation maps, and direct images.
 - The (un-blotted) segmentation map, used to derive the contamination models ( _e.g._ "_nis-wfss-ir\_seg.fits_" ).
 - The direct images used to create the reference catalogue ( _e.g._ "_nis-wfss-ir\_drz\_sci.fits_", "_nis-wfss-{f}\_drz\_sci.fits_" for filters _f_).

## Example

First, begin by importing the main class.

```python
from pygrife import GrismExtractor
run_app()
```
To instantiate the `GrismExtractor` object, we need to supply the `root_name` of the processed files, their current location (see [Necessary Files](#necessary-files)), and the output directory. `PyGriFE` copies the relevant files to the output directory, to preserve the originals and prevent any unexpected modifications.

```python
ge = GrismExtractor(
    field_root="nis-wfss",
    in_dir="probably/called/Prep",
    out_dir="use/a/new/name/for/this",
)
```
Next, we need to load the original segmentation map.

```python
ge.load_orig_seg_map(
    "/path/to/nis-wfss-ir_seg.fits",
)
```
We can now begin modifying the segmentation map to our liking. `ge.seg_map` is a 2D `numpy` array, and can be modified using any relevant operation. Several convenience methods are also built into `GrismExtractor`, such as the ability to read in a `regions` file written by [DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9).

```python
new_id = ge.add_reg_obj(
    reg_path="/some/path/to/targets_ds9_A2744_03.reg",
)
```
Alternatively, we could extract an object within a particular radius of a point (with many different ways to express this; see the documentation for further details).
```python
import numpy as np
import astropy.units as u

new_ids = np.append(
    new_id,
    ge.add_circ_obj(
        ra=3.61066,
        dec=-30.39560,
        unit="deg",
        radius=3*u.arcsec,
    ),
)
```
Now we've modified the segmentation map, we need to regenerate the object catalogue, and load the processed grism data.
```python
ge.regen_multiband_catalogue()
ge.load_grism_files(cpu_count=6)
```
Any objects can then be extracted, by supplying an array of IDs, and specifying the redshift range to search through.
```python
ge.extract_spectra(
    new_ids,
    z_range=[0.25,0.35],
)
```
If the sources are in a heavily contaminated part of the field, the contamination model can be refined using the best-fit spectra, and the process repeated.
```python
ge.refine_contam_model_with_fits(
    fit_files=[ge.out_dir / f"{ge.field_root}_{o_id:0>5}.full.fits" for o_id in new_ids],
)
ge.extract_spectra(
    new_ids,
    z_range=[0.25,0.35],
)
```
