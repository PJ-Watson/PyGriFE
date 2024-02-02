Installation
============

In all cases, it is strongly recommended to install ``PyGriFE`` into a new virtual environment, to minimise dependency conflicts (see [Requirements](#requirements)).

Building from Source
--------------------

To clone the latest GitHub repository, use this command:

.. code-block::

    git clone https://github.com/PJ-Watson/PyGriFE.git

To build and install ``PyGriFE``, run (from the root of the source tree):

.. code-block::

    python -m pip install .

Alternatively, if `grizli <https://github.com/gbrammer/grizli>`_ and `sep <https://github.com/kbarbary/sep>`_ are not already installed in the environment, they can be installed by running:

.. code-block::

    python -m pip install .[grizli,sep]

Necessary Files
---------------

``PyGriFE`` is designed to work with pre-processed grism data, specifically the output from `grizli <https://github.com/gbrammer/grizli>`_.
If you do not have the following files, ``PyGriFe`` will not be terribly useful:

* "*\*GrismFLT.fits*" and "*\*GrismFLT.pkl*" files, containing the grism exposures, contamination models, blotted segmentation maps, and direct images.
* The (un-blotted) segmentation map, used to derive the contamination models ( *e.g.* "*nis-wfss-ir\_seg.fits*" ).
* The direct images used to create the reference catalogue ( *e.g.* "*nis-wfss-ir\_drz\_sci.fits*", "*nis-wfss-{f}\_drz\_sci.fits*" for filters *f* ).

Requirements
------------

``PyGriFE`` has been tested with Python 3.10, and is developed primarily on Python 3.11.
The following packages are required for the most basic installation, and will be installed automatically if using ``pip``.

* `Python <https://www.python.org/>`_ 3.10 or later

* `Astropy <https://www.astropy.org/>`_ 5.3 or later

* `NumPy <https://www.numpy.org/>`_ 1.24 or later

* `Matplotlib <https://matplotlib.org/>`_ 3.6 or later

* `tqdm <https://tqdm.github.io/>`_ 4.66 or later

* `astropy-regions <https://astropy-regions.readthedocs.io>`_ 0.8 or later


Additional Dependencies
-----------------------

In addition to the packages listed above, ``PyGriFE`` requires several other packages to function correctly.
These are not installed by default to avoid conflicts with existing installations, but can be given as an option
to ``pip`` (see :ref:`installation:building from source`).

Grizli
^^^^^^
To extract objects, a working installation of `grizli <https://github.com/gbrammer/grizli>`_ is required.
``PyGriFE`` modifies some of the ``grizli`` class methods in order to extract arbitrary regions, and cannot
be guaranteed to work with every version (if ``PyGriFE`` encounters any compatibility problems, please raise
an issue `here <https://github.com/PJ-Watson/PyGriFE/issues>`_, rather than bothering the ``grizli`` developers).
The current tested version is ``grizli==1.11``.

SEP
^^^
Both ``PyGriFE`` and ``grizli`` rely on `SEP <https://github.com/kbarbary/sep>`_, the Python implementation of
`Source Extractor <http://www.astromatic.net/software/sextractor>`_ (`Bertin & Arnouts 1996 <http://adsabs.harvard.edu/abs/1996A%26AS..117..393B>`_).
Unfortunately, since the original repository no longer appears to be maintained, it is necessary to install the
fork maintained at `PJ-Watson/sep <https://github.com/PJ-Watson/sep>`_. This includes several bug fixes, and the
functionality necessary to rebuild a catalogue from an existing segmentation map.

GUI
^^^
To run the GUI, the following packages are required:
 - `PyQt6 <https://www.riverbankcomputing.com/software/pyqt/>`_ 6.6 or later
 - `qimage2ndarray <https://github.com/hmeine/qimage2ndarray>`_ 1.10 or later

To install the latest versions automatically during the build process for ``PyGriFE``, run

.. code-block::

    python -m pip install .[GUI]
