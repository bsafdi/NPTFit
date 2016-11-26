NPTFit Documentation
====================

NPTFit is a specialized Python/Cython package that implements
Non-Poissonian Template Fitting (NPTF), originally developed for characterizing 
populations of unresolved point sources. The main features of the package are

- Fast evaluation of likelihoods for NPTF analyses
- Easy-to-use interface for performing non-Poissonian (as well as standard Poissonian)
  template fits using MultiNest or other inference tools
- Ability to include an arbitrary number of point source
  templates, with an arbitrary number of degrees of freedom in the modeled flux distribution
- Modules for analyzing and plotting the results of an
  NPTF

The most up-to-date version of the code can be found at https://github.com/bsafdi/NPTFit.

Installation
------------

Out of the box, NPTFit relies on `MultiNest <https://ccpforge.cse.rl.ac.uk/gf/project/multinest/>`_ for Bayesian inference, which must be
installed and linked prior to use. 

NPTFit supports both Python 2 and 3, specifically 2.7 and 3.5. It may work with earlier 3.* versions, although this has not been tested.

Make sure Cython is installed (*e.g.* :code:`pip install Cython`). NPTFit along with it's dependent Python packages can then be installed with

.. code:: sh

  $ python setup.py install

which also builds the Cython modules. To just compile the Cython modules locally

.. code:: sh

  $ make build

The code is parallelizable through MPI (*e.g.* `OpenMPI <https://www.open-mpi.org/software/ompi/v2.0/>`_), which can
considerably speed up computationally intensive scans. This requires the MPI4py Python package for use with MultiNest, which
can be installed, for example, with ``pip``:

.. code:: sh

  $ pip install mpi4py

Contents
--------

.. toctree::
   :maxdepth: 2

   initializing_scan.rst

   loading_data_exposure.rst

   creating_masks.rst

   adding_models.rst

   PSFCorrection.rst

   running_the_scan.rst

   analyzing_results.rst

Example notebooks
-----------------

.. toctree::
   :maxdepth: 1

   Example1_Overview_of_the_Fermi_Data.rst

   Example2_Creating_Masks.rst

   Example3_Running_Poissonian_Scans.rst
   
   Example4_PSF_Correction.rst

   Example5_Running_nonPoissonian_Scans.rst

   Example6_Manual_nonPoissonian_Likelihood.rst

   Example7_Galactic_Center_nonPoissonian.rst

   Example8_Analysis.rst

   Example9_HighLat_Analysis.rst

.. TIP::
   When using NPTFit we recommend the use of National Pipe Thread (NPT) rather than the more common British Standard Pipe (BSP). For details, see `NPT <https://en.wikipedia.org/wiki/National_pipe_thread>`_. 

.. _``MultiNest``: https://ccpforge.cse.rl.ac.uk/gf/project/multinest/
.. _``OpenMPI``: https://www.open-mpi.org/software/ompi/v2.0/
.. _``NPT``: https://en.wikipedia.org/wiki/National_pipe_thread
