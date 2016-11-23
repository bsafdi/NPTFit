PSF correction
--------------

Fundamentally the presence of a non-zero PSF implies that the photons from any point source will be smeared out into some region around its true location. This effect must be accounted for in the NPTF. This is achieved via a function :math:`\rho(f)`. In the code we discretize :math:`\rho(f)` as an approximation to the full function.

If the angular reconstruction of the data is perfect, then :math:`\rho(f) = \delta(f-1)`. In many situations, such as for the Fermi data at higher energies, a Gaussian approximation of the PSF will suffice. Even then there are a number of variables that go into evaluating the correction, as shown below. Finally we will show how the code can be used for the case of non-Gaussian PSFs.

The module ``psf_correction`` is used to account for the instrument point spread function (PSF). 

By default, the correction for a Gaussian PSF can be determined as follows:

.. code:: python

    >>> import psf_correction as pc
    >>> pc_inst = pc.psf_correction(psf_sigma_deg)
    >>> f_ary, df_rho_div_f_ary = pc_inst.f_ary, pc_inst.df_rho_div_f_ary

where ``psf_sigma_deg`` is the 68% containment of the Gaussian PSF.

The two outputs of an instance of psf_correction are: ``f_ary``, an array of :math:`f` values; and ``df_rho_div_f_ary``, an associated array of :math:`\Delta f \rho(f)/f` values, where :math:`\Delta f` is the width of the ``f_ary`` bins.

Additional PSF options:

+---------------------+-----------+---------------------------------------------------------+
| Argument            | Defaults  | Purpose                                                 |
+=====================+===========+=========================================================+
| ``psf_dir``         | ``None``  | Directory where PSF files are stored                    |
+---------------------+-----------+---------------------------------------------------------+
| ``num_f_bins``      | 10        | Number of bins to calculate rho(f) in                   |
+---------------------+-----------+---------------------------------------------------------+
| ``n_psf``           | 50000     | Number of psfs to place when calculating rho(f)         |
+---------------------+-----------+---------------------------------------------------------+
| ``n_pts_per_psf``   | 1000      | Number of points to place per PSF in calculation        |
+---------------------+-----------+---------------------------------------------------------+
| ``f_trunc``         | 0.01      | Minimum flux fraction to keep track of                  |
+---------------------+-----------+---------------------------------------------------------+
| ``psf_sigma_deg``   | ``None``  | 68% containment of the Gaussian PSF                     |
+---------------------+-----------+---------------------------------------------------------+
| ``delay_compute``   | ``False`` | Set to true to define a custom PSF                      |
+---------------------+-----------+---------------------------------------------------------+
| ``nside``           | 128       | Pixelization of the data map                            |
+---------------------+-----------+---------------------------------------------------------+

We show how to pass the PSF information into NPTFit in :doc:`running_the_scan`. 

Using a custom PSF
~~~~~~~~~~~~~~~~~~

To implement a custom (non-Gaussian) PSF create an instance of
``psf_correction`` with the above parameters (ignore ``psf_sigma_deg``)
setting ``delay_compute`` = True.

Then redefine the following attributes:

+--------------------+------------------+-----------------------------------+
| Argument           | Defaults         | Purpose                           |
+====================+==================+===================================+
| ``psf_r_func``     | Gaussian         | The PSF as a function of r,       |
|                    |                  | distance from the center          |
+--------------------+------------------+-----------------------------------+
| ``sample_psf_max`` | ``5.*psf_sigma`` | Maximum distance to sample the    |
|                    |                  | PSF from the center,              | 
|                    |                  |                                   |
|                    |                  | should be larger for PSFs with    |
|                    |                  | long tails                        |
+--------------------+------------------+-----------------------------------+
| ``psf_samples``    | 10000            | Number of samples to make from    |
|                    |                  | the PSF (linearly                 |
|                    |                  |                                   |
|                    |                  | spaced) from 0 to                 |
|                    |                  | ``sample_psf_max``, should be     |
|                    |                  |                                   |
|                    |                  | large enough to adequately        |
|                    |                  | represent the full PSF            |
+--------------------+------------------+-----------------------------------+
| ``psf_tag``        | \-               | Label the PSF is saved with       |
+--------------------+------------------+-----------------------------------+


Finally, execute the calculation using

.. code:: python

    >>> pc_inst.make_or_load_psf_corr()

and the output can be extracted as in the Gaussian case.

.. NOTE::
   As the calculation of :math:`\rho(f)` can be time consuming, we always save the output to avoid recomputing the same correction twice.

.. WARNING::
   The PSF correction assumes it will be applied to a 
   `HEALPix <http://healpix.jpl.nasa.gov/>`_
   map of binning ``nside``. 

.. TIP::
   See :doc:`Example4_PSF_Correction` for an exposition of the options described here.

.. _``HEALPix``: http://healpix.jpl.nasa.gov/
