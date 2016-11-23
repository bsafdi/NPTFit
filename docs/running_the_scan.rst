Running the scan
----------------

The scan can be configured as follows:

.. code:: python

    >>> nptf.configure_for_scan(f_ary, df_rho_div_f_ary, nexp)
    
+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| Argument                | Default| Purpose                                                                                          |
+=========================+========+==================================================================================================+
| ``f_ary``               | ``[1]``| See :doc:`psf_correction`                                                                        |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| ``df_rho_div_f_ary``    | ``[1]``| See :doc:`psf_correction`                                                                        |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| ``nexp``                | 1      |  Number of exposure regions the scan is divided into when performing                             |
|                         |        |  into when performing an NPTF                                                                    |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+

.. NOTE::
   For no PSF correction, simply leave ``f_ary`` and ``df_rho_div_f_ary`` to their default values.

.. NOTE::
   If no non-Poissonian models are added when ``configure_for_scan`` is called, the
   likelihood will default to a pure Poissonian scan.

From here the scan can be run using:

.. code:: python

    >>> nptf.perform_scan(run_tag='example_run', nlive, pymultinest_options=None)

Scan options:

+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| Argument                | Default| Purpose                                                                                          |
+=========================+========+==================================================================================================+
| ``run_tag``             | None   | An optional custom label for the folder within ``work_dir/chains/tag``                           |
|                         |        |                                                                                                  |
|                         |        | where the MultiNest output is stored                                                             |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| ``nlive``               | 100    | Number of live sampling points in MultiNest                                                      |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+
| ``pymultinest_options`` | None   |  A custom set of MCMC options for ``pymultinest`` (see below)                                    |
+-------------------------+--------+--------------------------------------------------------------------------------------------------+

Setting custom MultiNest options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A custom set of MCMC options for ``pymultinest`` can be provided to
``nptfit.NPTF.perform_scan()`` with the argument ``pymultinest_options``. This
should be a dictionary of the form ``{'option': value ...}``. The
default options are:

.. code:: python

    pymultinest_options = {importance_nested_sampling: False, 
                           resume: False, verbose: True, 
                           make_or_load_psf_corr()sampling_efficiency: model, 
                           init_MPI: False, evidence_tolerance: 0.5,
                           const_efficiency_mode: False}

Parallel implementation
~~~~~~~~~~~~~~~~~~~~~~~

MultiNest uses MPI for parallel sampling, which can significantly speed up the NPTF. This can 
be seamlessly implemented as with ``NPTFit``. For example, to run on 12 CPU cores:

::

    >>> mpirun -np 12 nptf_run.py
