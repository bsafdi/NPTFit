
Example 8: Analyzing the Results of an NPTFit Run
=================================================

While the chain samples of a non-Poissonian fit performed using
MultiNest can be readily accessed, we provide a basic analysis module,
``dnds_analysis.py`` that contains helper functions to: 1. Make triangle
plots; 2. Get template intensities; 3. Plot source count distributions;
4. Plot flux fractions; 5. Access individual posteriors; and 6. Get
Bayesian log-evidences.

In this example we provide the details of how to use each function.

**NB:** Example 7 must be run before this notebook. Note that the run
performed there was with a low nside and a fixed diffuse model, so the
results below should be interpreted only as approximate.

.. code:: python

    # Import relevant modules
    
    %matplotlib inline
    %load_ext autoreload
    %autoreload 2
    
    import numpy as np
    import corner
    import matplotlib.pyplot as plt
    
    from NPTFit import nptfit # module for performing scan
    from NPTFit import create_mask as cm # module for creating the mask
    from NPTFit import dnds_analysis # module for analysing the output
    from NPTFit import psf_correction as pc # module for determining the PSF correction
    
    from __future__ import print_function

Analysis
--------

At the outset, an instance of ``nptfit.NPTF`` must be created as was
done when initiating and performing the scan. The process here is the
same as in Example 7, up to configuring the scan. Finally, the scan is
loaded with ``n.load_scan()``.

.. code:: python

    n = nptfit.NPTF(tag='GCE_Example_newf')

.. code:: python

    fermi_data = np.load('fermi_data/fermidata_counts.npy')
    fermi_exposure = np.load('fermi_data/fermidata_exposure.npy')
    n.load_data(fermi_data, fermi_exposure)

.. code:: python

    pscmask=np.array(np.load('fermi_data/fermi_pscmask_cons.npy'), dtype=bool)
    analysis_mask = cm.make_mask_total(band_mask = True, band_mask_range = 2,
                                       mask_ring = True, inner = 0, outer = 30,
                                       custom_mask = pscmask)
    n.load_mask(analysis_mask)

.. code:: python

    dif = np.load('fermi_data/template_dif.npy')
    iso = np.load('fermi_data/template_iso.npy')
    bub = np.load('fermi_data/template_bub.npy')
    gce = np.load('fermi_data/template_gce.npy')
    dsk = np.load('fermi_data/template_dsk.npy')
    
    n.add_template(dif, 'dif')
    n.add_template(iso, 'iso')
    n.add_template(bub, 'bub')
    n.add_template(gce, 'gce')
    n.add_template(dsk, 'dsk')

.. code:: python

    n.add_poiss_model('dif', '$A_\mathrm{dif}$', fixed=True, fixed_norm=14.67)
    n.add_poiss_model('iso', '$A_\mathrm{iso}$', [0,2], False)
    n.add_poiss_model('gce', '$A_\mathrm{gce}$', [0,2], False)
    n.add_poiss_model('bub', '$A_\mathrm{bub}$', [0,2], False)

.. code:: python

    n.add_non_poiss_model('gce',
                          ['$A_\mathrm{gce}^\mathrm{ps}$','$n_1^\mathrm{gce}$','$n_2^\mathrm{gce}$','$S_b^{(1), \mathrm{gce}}$'],
                          [[-6,1],[2.05,30],[-2,1.95],[0.05,40]],
                          [True,False,False,False])
    n.add_non_poiss_model('dsk',
                          ['$A_\mathrm{dsk}^\mathrm{ps}$','$n_1^\mathrm{dsk}$','$n_2^\mathrm{dsk}$','$S_b^{(1), \mathrm{dsk}}$'],
                          [[-6,1],[2.05,30],[-2,1.95],[0.05,40]],
                          [True,False,False,False])

.. code:: python

    pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812)
    f_ary, df_rho_div_f_ary = pc_inst.f_ary, pc_inst.df_rho_div_f_ary


.. parsed-literal::

    Loading the psf correction from: /group/hepheno/smsharma/NPTFit/examples/psf_dir/gauss_128_0.181_10_50000_1000_0.01.npy


.. code:: python

    n.configure_for_scan(f_ary, df_rho_div_f_ary, nexp=1)


.. parsed-literal::

    The number of parameters to be fit is 11


Finally, instead of running the scan we simply load the completed scan
performed in Example 7.

.. code:: python

    n.load_scan()


.. parsed-literal::

      analysing data from /group/hepheno/smsharma/NPTFit/examples/chains/GCE_Example_newf/.txt


Analysis
--------

An instance of ``nptf.NPTF`` with a loaded scan as above can already be
used to access the posterior chains with ``n.samples``:

.. code:: python

    print(np.shape(n.samples))
    print(n.samples)


.. parsed-literal::

    (741, 11)
    [[  3.64691822e-01   1.34740973e-01   9.31172758e-01 ...,   2.20491850e+01
        5.44277104e-01   3.02346663e+01]
     [  6.43735503e-01   4.23342456e-02   7.53877709e-01 ...,   2.09548864e+01
       -7.15130115e-02   3.08319347e+01]
     [  4.26872933e-01   9.17816072e-02   9.68932781e-01 ...,   1.92173721e+01
       -1.87495053e+00   2.68394950e+01]
     ..., 
     [  4.43485759e-01   2.29078417e-03   8.81020387e-01 ...,   1.28963036e+01
       -1.84868559e+00   3.08482690e+01]
     [  4.55746431e-01   1.41377732e-02   9.16967489e-01 ...,   2.78378738e+01
       -1.39216537e+00   2.98667420e+01]
     [  4.81895556e-01   1.80832165e-02   8.79153360e-01 ...,   2.43908162e+01
       -5.38391016e-01   3.00742214e+01]]


In the analysis module described next we provide basic helper functions
to load in and manipulate these chain samples.

0. Initialize Analysis Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first thing to do is initialize an instance of the analysis module,
``dnds_analysis`` from ``dnds_analysis.py`` with a provided instance of
``nptfit.NPTF``. The ``NPTF`` instance should have a scan already loaded
in, as done with ``n.load_scan()`` above.

.. code:: python

    an = dnds_analysis.Analysis(n)

``dnds_analysis`` has an optional argument ``mask``, which if unset
defaults to the mask in the passed instance of ``NPTF``. If a mask is
given, however, then the analysis will be performed in a different ROI
to the main run.

1. Make triangle plots
~~~~~~~~~~~~~~~~~~~~~~

Triangle/corner plots let us visualize multidimensional samples using a
scatterplot matrix. A triangle plot with the default options can be made
as follows.

.. code:: python

    an.make_triangle()



.. image:: Example8_Analysis_files/Example8_Analysis_25_0.png


To use your own custom plotting options, use corner as follows

.. code:: python

    corner.corner(an.nptf.samples, labels=an.nptf.params, range=[1 for i in range(an.nptf.n_params)])

with additional arguments as specified in
http://corner.readthedocs.io/en/latest/.

2. Get Intensities
~~~~~~~~~~~~~~~~~~

Template intensities can be calculated with

.. code:: python

    dnds_analysis.return_intensity_arrays_poiss(comp)
    dnds_analysis.return_intensity_arrays_non_poiss(comp)

for the Poissonian and non-Poissonian templates respectively. This
returns an intensity array corresponding to each chain sample associated
with the template ``comp``.

The NPT intensity is calculated by integrating up
:math:`\int_{S_{min}}^{S_{max}} dS~S~dN/dS`. This is approximated as a
sum between :math:`S_{min}` and :math:`S_{max}`. The options associated
with the non-Poissonian template intensity are:

+--------------+-----------------+--------------------------------------------+
| Argument     | Default Value   | Purpose                                    |
+==============+=================+============================================+
| ``comp``     | -               | The NPT key                                |
+--------------+-----------------+--------------------------------------------+
| ``smin``     | 0.01            | Minimum counts to sum up from              |
+--------------+-----------------+--------------------------------------------+
| ``smax``     | 10000           | Maximum counts to sum up to                |
+--------------+-----------------+--------------------------------------------+
| ``nsteps``   | 10000           | Number of bins in ``s`` while summing up   |
+--------------+-----------------+--------------------------------------------+

We can then look at the quantiles of this distribution, for example to
see the middle 68% along with the medians of the GCE and disk NPT as
well as that of the GCE PT:

.. code:: python

    print("GCE NPT Intensity", corner.quantile(an.return_intensity_arrays_non_poiss('gce'),[0.16,0.5,0.84]), "ph/cm^2/s")
    print("Disk NPT Intensity", corner.quantile(an.return_intensity_arrays_non_poiss('dsk'),[0.16,0.5,0.84]), "ph/cm^2/s")
    print("GCE PT Intensity", corner.quantile(an.return_intensity_arrays_poiss('gce'),[0.16,0.5,0.84]), "ph/cm^2/s")


.. parsed-literal::

    GCE NPT Intensity [  4.84629246e-08   6.40066569e-08   7.72576661e-08] ph/cm^2/s
    Disk NPT Intensity [  4.63451387e-08   6.27115423e-08   8.07192577e-08] ph/cm^2/s
    GCE PT Intensity [  8.98465954e-10   3.19275586e-09   8.24976317e-09] ph/cm^2/s


3. Plot Source Count Distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The posterior arrays for the source count distributions :math:`dN/dF`
[counts:math:`^{-1}` cm\ :math:`^2` s deg\ :math:`^{-2}`] associated
with a given template ``comp`` at a given ``flux``
[counts/cm:math:`^2`/s] can be obtained using

.. code:: python

    dnds.return_dndf_arrays(comp,flux)

The quantiles of this can then be obtained as before. For example, the
middle 68% and medians for the GCE and disk non-Poissonian templates:

.. code:: python

    print(corner.quantile(an.return_dndf_arrays('gce',1e-12),[0.16,0.5,0.84]))
    print(corner.quantile(an.return_dndf_arrays('dsk',1e-12),[0.16,0.5,0.84]))


.. parsed-literal::

    [  1.04149898e+05   3.60411984e+06   1.71299721e+08]
    [  3.51206262e+04   1.32927373e+06   1.28748409e+08]


The following arrays are used to show the resolved 3FGL points sources
and associated Poisson errors as appropriate for the plots below. For
how these were obtained, see `this
snippet <https://gist.github.com/smsharma/829296c483a92528ab8bbba0d1439e88>`__.

.. code:: python

    x_counts, y_counts, error_L, error_H, x_errors_L, x_errors_H = \
    [np.array([  1.36887451e-10,   2.56502091e-10,   4.80638086e-10,
              9.00628020e-10,   1.68761248e-09,   3.16227766e-09,
              5.92553098e-09,   1.11033632e-08,   2.08056754e-08,
              3.89860370e-08,   7.30527154e-08]),
     np.array([  1.04000127e+08,   1.83397053e+08,   9.65856820e+07,
              1.51198295e+07,   4.76804443e+06,   9.78677656e+05,
              2.08916332e+05,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00,   0.00000000e+00]),
     np.array([  2.14237668e+07,   2.08831658e+07,   1.10708578e+07,
              3.18362798e+06,   1.29929969e+06,   4.21069315e+05,
              1.34538182e+05,  -5.57461814e-04,  -2.97500603e-04,
             -1.58767124e-04,  -8.47292389e-05]),
     np.array([  2.63822671e+07,   2.34164673e+07,   1.24232945e+07,
              3.93887993e+06,   1.71404939e+06,   6.58746511e+05,
              2.74201404e+05,   1.02159419e+05,   5.45194091e+04,
              2.90953689e+04,   1.55273233e+04]),
     np.array([  3.68874510e-11,   6.91203483e-11,   1.29518913e-10,
              2.42694796e-10,   4.54765736e-10,   8.52147960e-10,
              1.59676969e-09,   2.99205487e-09,   5.60656455e-09,
              1.05056783e-08,   1.96857231e-08]),
     np.array([  5.04942913e-11,   9.46170829e-11,   1.77295138e-10,
              3.32218719e-10,   6.22517224e-10,   1.16648362e-09,
              2.18577733e-09,   4.09574765e-09,   7.67468330e-09,
              1.43809553e-08,   2.69472846e-08])]

The source count distribution can be plotted with

.. code:: python

    dnds.plot_source_count_median(comp, smin, smax, nsteps, spow, **kwargs)
    dnds.plot_source_count_band(comp, smin, smax, nsteps, spow, qs, **kwargs)

The options being the same as for obtaining the NPT intensity above.
Additionally, ``spow`` is the power :math:`n` in :math:`F^ndN/dF` to
return while plotting, and ``qs`` is an array of quantiles for which to
return the dN/dF band. We plot here the median in addition to 68% and
95% confidence intervals.

.. code:: python

    plt.figure(figsize=[6,5])
    
    an.plot_source_count_median('dsk',smin=0.01,smax=1000,nsteps=1000,color='royalblue',spow=2,label='Disk PS')
    an.plot_source_count_band('dsk',smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],color='royalblue',alpha=0.15,spow=2)
    an.plot_source_count_band('dsk',smin=0.01,smax=1000,nsteps=1000,qs=[0.025,0.5,0.975],color='royalblue',alpha=0.1,spow=2)
    
    
    an.plot_source_count_median('gce',smin=0.01,smax=1000,nsteps=1000,color='firebrick',spow=2,label='GCE PS')
    an.plot_source_count_band('gce',smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],color='firebrick',alpha=0.15,spow=2)
    an.plot_source_count_band('gce',smin=0.01,smax=1000,nsteps=1000,qs=[0.025,0.5,0.975],color='firebrick',alpha=0.1,spow=2)
    
    plt.errorbar(x_counts,x_counts**2*y_counts,xerr=[x_errors_L,x_errors_H],yerr=x_counts**2*np.array([error_L,error_H]), fmt='o', color='black', label='3FGL PS')
    
    
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([1e-10,1e-8])
    plt.ylim([2e-13,1e-10])
    
    plt.tick_params(axis='x', length=5, width=2, labelsize=18)
    plt.tick_params(axis='y', length=5, width=2, labelsize=18)
    plt.ylabel('$F^2 dN/dF$ [counts cm$^{-2}$s$^{-1}$deg$^{-2}$]', fontsize=18)
    plt.xlabel('$F$  [counts cm$^{-2}$ s$^{-1}$]', fontsize=18)
    plt.title(r'Galactic Center NPTF', y=1.02)
    plt.legend(fancybox=True, loc='lower right');
    plt.tight_layout()
    
    
    # plt.savefig("dnds_masked.pdf")



.. image:: Example8_Analysis_files/Example8_Analysis_36_0.png


As some references also show :math:`dN/dF`, and we give an example of it
below, also demonstrating the use of ``spow``.

.. code:: python

    plt.figure(figsize=[6,5])
    
    an.plot_source_count_median('dsk',smin=0.01,smax=1000,nsteps=1000,color='royalblue',spow=0,label='Disk PS')
    an.plot_source_count_band('dsk',smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],color='royalblue',alpha=0.15,spow=0)
    an.plot_source_count_band('dsk',smin=0.01,smax=1000,nsteps=1000,qs=[0.025,0.5,0.975],color='royalblue',alpha=0.1,spow=0)
    
    
    an.plot_source_count_median('gce',smin=0.01,smax=1000,nsteps=1000,color='firebrick',spow=0,label='GCE PS')
    an.plot_source_count_band('gce',smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],color='firebrick',alpha=0.15,spow=0)
    an.plot_source_count_band('gce',smin=0.01,smax=1000,nsteps=1000,qs=[0.025,0.5,0.975],color='firebrick',alpha=0.1,spow=0)
    
    plt.errorbar(x_counts, y_counts,xerr=[x_errors_L,x_errors_H],yerr=np.array([error_L,error_H]), fmt='o', color='black', label='3FGL PS')
    
    
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([5e-11,5e-9])
    plt.ylim([2e5,2e9])
    plt.tick_params(axis='x', length=5, width=2, labelsize=18)
    plt.tick_params(axis='y', length=5, width=2, labelsize=18)
    plt.ylabel('$dN/dF$ [counts$^{-1}$cm$^2$ s deg$^{-2}$]', fontsize=18)
    plt.xlabel('$F$  [counts cm$^{-2}$ s$^{-1}$]', fontsize=18)
    plt.title('Galactic Center NPTF', y=1.02)
    plt.legend(fancybox=True);



.. image:: Example8_Analysis_files/Example8_Analysis_38_0.png


4. Plot Intensity Fractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intensity fractions (fraction of template intensity to total intensity)
for Poissonian and non-Poissonian templates respectively can be plotting
using

.. code:: python

    dnds.plot_intensity_fraction_poiss(comp, bins, **kwargs)
    dnds.plot_intensity_fraction_non_poiss(comp, bins, **kwargs)

where ``comp`` is the template key, ``bins`` is the number of bins
between 0 and 100 and ``**kwargs`` specify plotting options.

.. code:: python

    an.plot_intensity_fraction_non_poiss('gce', bins=800, color='firebrick', label='GCE PS')
    an.plot_intensity_fraction_poiss('gce', bins=800, color='forestgreen', label='GCE DM')
    plt.xlabel('Flux fraction (%)')
    plt.legend(fancybox = True)
    plt.xlim(0,6);
    plt.ylim(0,.3);



.. image:: Example8_Analysis_files/Example8_Analysis_41_0.png


This plot makes it clear, that when given the choice, the fit prefers to
put the GCE flux into point sources rather than diffuse emission.

5. Access Parameter Posteriors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the posteriors can be accessed with ``n.samples`` (or
``an.nptf.samples``) as above, the following functions provide a useful
interfact to access individual parameters:

.. code:: python

    dnds_analysis.return_poiss_parameter_posteriors(comp)
    dnds_analysis.return_poiss_parameter_posteriors(comp)

where ``comp`` is the (non-)Poissonian template key.

Poissonian parameters
^^^^^^^^^^^^^^^^^^^^^

Posterior normalizations of Poissonian parameters can be loaded simply
as:

.. code:: python

    Aiso_poiss_post = an.return_poiss_parameter_posteriors('iso')
    Agce_poiss_post = an.return_poiss_parameter_posteriors('gce')
    Abub_poiss_post = an.return_poiss_parameter_posteriors('bub')

These can then be use in any way required, for example simply plotted:

.. code:: python

    f, axarr = plt.subplots(nrows = 1, ncols=3)
    f.set_figwidth(12)
    f.set_figheight(4)
    
    axarr[0].hist(Aiso_poiss_post, histtype='stepfilled', color='cornflowerblue', bins=np.linspace(0,1.,30), alpha=.4);
    axarr[0].set_title('$A_\mathrm{iso}$')
    axarr[1].hist(Agce_poiss_post, histtype='stepfilled', color='lightsalmon', bins=np.linspace(0,.2,30), alpha=.4);
    axarr[1].set_title('$A_\mathrm{gce}$')
    axarr[2].hist(Abub_poiss_post, histtype='stepfilled', color='plum', bins=np.linspace(.5,1.5,30), alpha=.4);
    axarr[2].set_title('$A_\mathrm{bub}$')
    
    plt.setp([a.get_yticklabels() for a in axarr], visible=False);
    
    plt.tight_layout()



.. image:: Example8_Analysis_files/Example8_Analysis_49_0.png


Non-poissonian parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

A similar syntax can be used to extract the non-Poissonian parameters.

.. code:: python

    Agce_non_poiss_post, n_non_poiss_post, Sb_non_poiss_post = an.return_non_poiss_parameter_posteriors('gce')

.. code:: python

    f, axarr = plt.subplots(2, 2);
    f.set_figwidth(8)
    f.set_figheight(8)
    
    
    axarr[0, 0].hist(Agce_non_poiss_post, histtype='stepfilled', color='cornflowerblue', bins=np.linspace(0,0.02,30), alpha=.4);
    axarr[0, 0].set_title('$A_\mathrm{gce}^\mathrm{ps}$')
    axarr[0, 1].hist(n_non_poiss_post[0], histtype='stepfilled', color='lightsalmon', bins=np.linspace(2,30,30), alpha=.4);
    axarr[0, 1].set_title('$n_1^\mathrm{gce}$')
    axarr[1, 0].hist(n_non_poiss_post[1], histtype='stepfilled', color='lightsalmon', bins=np.linspace(-2,2,30), alpha=.4);
    axarr[1, 0].set_title('$n_2^\mathrm{gce}$')
    axarr[1, 1].hist(Sb_non_poiss_post, histtype='stepfilled', color='plum', bins=np.linspace(0,40,30), alpha=.4);
    axarr[1, 1].set_title('$S_b^{(1), \mathrm{gce}}$')
    
    plt.setp(axarr[0, 0], xticks=[x*0.01 for x in range(5)])
    plt.setp(axarr[1, 0], xticks=[x*1.0-2 for x in range(5)])
    plt.setp(axarr[1, 1], xticks=[x*10 for x in range(6)])
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False);
    plt.setp([a.get_yticklabels() for a in axarr[:, 0]], visible=False);
    
    plt.tight_layout()



.. image:: Example8_Analysis_files/Example8_Analysis_53_0.png


6. Bayesian log-evidence
~~~~~~~~~~~~~~~~~~~~~~~~

Finally the Bayesian log-evidence and associated error can be accessed
as follows.

.. code:: python

    l_be, l_be_err = an.get_log_evidence()
    print(l_be, l_be_err)


.. parsed-literal::

    -29466.494212 0.411995897483

