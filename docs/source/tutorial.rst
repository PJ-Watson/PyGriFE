Basic Usage
===========

There are two methods to run ``PyGriFE``. The recommended method is to use a text-based python session,
although there is an optional GUI.

Python
------

First, begin by importing the main class into the python environment.

.. code-block:: python

    from pygrife import GrismExtractor

To instantiate the `~pygrife.extractor_main.GrismExtractor` object, we need to supply the ``root_name`` of the processed files, their current location (see necessary files), and the output directory. `~pygrife.extractor_main.GrismExtractor` copies the relevant files to the output directory, to preserve the originals and prevent any unexpected modifications.

.. code-block:: python

    ge = GrismExtractor(
        field_root="nis-wfss",
        in_dir="probably/called/Prep",
        out_dir="use/a/new/name/for/this",
    )

Next, we need to load the original segmentation map.

.. code-block:: python

    ge.load_orig_seg_map(
        "/path/to/nis-wfss-ir_seg.fits",
    )

We can now begin modifying the segmentation map to our liking. ``ge.seg_map`` is a 2D `numpy` array, and can be modified using any relevant operation. Several convenience methods are also built into `~pygrife.extractor_main.GrismExtractor`, such as the ability to read in a "regions" file written by `DS9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_.

.. code-block:: python

    new_id = ge.add_reg_obj(
        reg_path="/some/path/to/targets_ds9_A2744_03.reg",
    )

Alternatively, we could extract an object within a particular radius of a point (with many different ways to express this; see the documentation of `~pygrife.extractor_main.GrismExtractor.add_circ_obj` for further details).

.. code-block:: python

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

Now we've modified the segmentation map, we need to regenerate the object catalogue, and load the processed grism data.

.. code-block:: python

    ge.regen_multiband_catalogue()
    ge.load_grism_files(cpu_count=6)

Any objects can then be extracted, by supplying an array of IDs, and specifying the redshift range to search through.

.. code-block:: python

    ge.extract_spectra(
        new_ids,
        z_range=[0.25,0.35],
    )

If you've made it this far without any problems, congratulations! You should now have all the beam cutouts,
extracted spectra, and other paraphernalia associated with your objects.
Of course, if the sources are in a heavily contaminated part of the field, the contamination model can be
refined using the best-fit spectra, and the process repeated until a satisfactory extraction is obtained.

.. code-block:: python

    ge.refine_contam_model_with_fits(
        fit_files=[ge.out_dir / f"{ge.field_root}_{o_id:0>5}.full.fits" for o_id in new_ids],
    )
    ge.extract_spectra(
        new_ids,
        z_range=[0.25,0.35],
    )

GUI
---

``PyGriFE`` also includes a graphical interface, to view and modify the segmentation map.
This feature is under active development, but can be run after installing the optional
dependencies (see :ref:`Installation: GUI<installation:gui>` for more details).
The GUI itself can be used by running the following.

.. code-block:: python

    from pygrife.GUI import run_GUI

    run_GUI()

Note that even when feature complete, at least one thread will be dedicated to running the GUI,
and so the method detailed above in :ref:`Basic Usage: Python<tutorial:python>` will be slightly more performant.
