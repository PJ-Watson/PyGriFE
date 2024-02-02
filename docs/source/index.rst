.. PyGriFE documentation master file, created by
   sphinx-quickstart on Thu Feb  1 13:00:22 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyGriFE
=======

A package to aid in the forced extraction of pre-processed grism data.

About
-----

This package was originally written for the `GLASS-JWST collaboration <https://glass.astro.ucla.edu/ers/>`_, but
should be more widely applicable to any slitless spectroscopic observations.
The intended purpose is for observations which have already been processed in some form,
and where the user wishes to extract a particular object or region, such as the very centre of a large galaxy, or
only the extended tail of a ram-pressure-stripped galaxy.
For an example of basic usage, see :doc:`tutorial`, or for the complete API documentation, see :doc:`modules`.

This package builds upon the output files and code used in the
`Grism redshift & line analysis software, grizli <https://grizli.readthedocs.io/en/latest/index.html>`_,
by `G. Brammer et al. <https://doi.org/10.5281/zenodo.1146904>`_.
Whilst it is not strictly necessary to master every technical detail of ``grizli`` in order to use
``PyGriFE``, a basic understanding of the process is invaluable.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   modules


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
