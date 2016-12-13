NPTFit
======

**Non-Poissonian template fitting in Python/Cython**

|Build Status| |Code Coverage| |Documentation Status| |License: MIT|

AUTHORS
-------

-  Benjamin Safdi; bsafdi at mit dot edu
-  Nicholas Rodd; nrodd at mit dot edu
-  Siddharth Mishra-Sharma; smsharma at princeton dot edu

A full list of the people who have contributed to NPTFit can be found
here:
`AUTHORS.txt <https://github.com/bsafdi/NPTFit/blob/master/AUTHORS.txt>`__.

If you make use of NPTFit in a publication, please cite
`1612.03173 <https://arxiv.org/abs/1612.03173>`__.

INSTALLATION AND DEPENDENCIES
-----------------------------


Out of the box, NPTFit relies on `MultiNest <https://ccpforge.cse.rl.ac.uk/gf/project/multinest/>`_ for Bayesian inference, which must be
installed and linked prior to use. See `here <http://monte-python.readthedocs.io/en/latest/nested.html>`_ for a helpful set of installation instructions.

NPTFit supports both Python 2 and 3, specifically 2.7 and 3.5. It may work with earlier 3.* versions, although this has not been tested.

Make sure Cython is installed (*e.g.* :code:`pip install Cython`). The easiest way to install NPTFit along with it's dependent Python packages 
is using ``pip``:

.. code:: sh

  $ pip install NPTFit

or using the setup script:

.. code:: sh

  $ python setup.py install

which also builds the Cython modules. To just compile the Cython modules locally:

.. code:: sh

  $ make build

The code is parallelizable through MPI (*e.g.* `OpenMPI <https://www.open-mpi.org/software/ompi/v2.0/>`_), which can
considerably speed up computationally intensive scans. This requires the MPI4Py Python package for use with MultiNest, which
can be installed, for example, with ``pip``:


.. code:: sh

  $ pip install mpi4py

DOCUMENTATION AND EXAMPLES
--------------------------

Detailed documentation of the code can be found
`here <http://nptfit.readthedocs.io/en/latest/>`__. There we also
provide a series of examples, all of which are also available either in
the form of interactive Jupyter notebooks, found
`here <https://github.com/bsafdi/NPTFit/tree/master/examples>`__.

BASIC USAGE
-----------

Here's how easy it is to perform a non-Poissonian template fit with
NPTFit.

.. code:: Python


    # Import modules
    import numpy as np
    from NPTFit import nptfit # module for performing scan
    from NPTFit import create_mask as cm # module for creating the mask
    from NPTFit import psf_correction as pc # module for determining the PSF correction
    from NPTFit import dnds_analysis # module for analysing the output

    # Initiate NPTF
    n = nptfit.NPTF()

    # Load data and templates
    fermi_data = np.load('fermi_data/fermidata_counts.npy')
    fermi_exposure = np.load('fermi_data/fermidata_exposure.npy')
    iso_temp = np.load('fermi_data/template_iso.npy')

    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso_temp, 'iso')

    # Define the Region of Interest with a Mask
    analysis_mask = cm.make_mask_total(mask_ring=True, inner=0, outer=5, ring_b=90, ring_l=0)
    n.load_mask(analysis_mask)

    # Add a Poissonian and non-Poissonian model
    n.add_poiss_model('iso','$A_\mathrm{iso}$', False, fixed=True, fixed_norm=1.47)
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_b$'],
                          [[-6,1],[2.05,30],[-2,1.95]],
                          [True,False,False],
                          fixed_params = [[3,22.]])

    # Calculate the PSF Correction
    pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812)
    f_ary = pc_inst.f_ary
    df_rho_div_f_ary = pc_inst.df_rho_div_f_ary

    # Configure and perform scan
    n.configure_for_scan(f_ary=f_ary, df_rho_div_f_ary=df_rho_div_f_ary)
    n.perform_scan(nlive=500)

An interactive version of this example can be found in the example
`here <https://github.com/bsafdi/NPTFit/blob/master/examples/Example5_Running_nonPoissonian_Scans.ipynb>`__.

The following source-count distribution is an unmasked version of the one produced in
`this <https://github.com/bsafdi/NPTFit/blob/master/examples/Example8_Analysis.ipynb>`__
example, which analyzed the output of
`this <https://github.com/bsafdi/NPTFit/blob/master/examples/Example7_Galactic_Center_nonPoissonian.ipynb>`__
example exploring the point source origin of the galactic center excess.

.. figure:: https://github.com/bsafdi/NPTFit/blob/master/docs/GCE_unmasked.png
   :alt: SourceCount

ISSUES
------

Problems with the code should be reported to the authors, or preferably
noted through the `issue
tracker <https://github.com/bsafdi/NPTFit/issues>`__.

.. |Code Coverage| image:: https://codecov.io/gh/bsafdi/NPTFit/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/bsafdi/NPTFit
.. |Build Status| image:: https://travis-ci.org/bsafdi/NPTFit.svg?branch=master
   :target: https://travis-ci.org/bsafdi/NPTFit
.. |Documentation Status| image:: https://readthedocs.org/projects/nptfit/badge/?version=latest
   :target: http://nptfit.readthedocs.io/en/latest/?badge=latest
.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT

