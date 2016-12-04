Analyzing results of a scan
---------------------------

.. TIP::
   See the :doc:`Example8_Analysis` for an interactive version of the analysis options described here.


While the chain samples of a non-Poissonian fit performed using MultiNest can be readily 
accessed, a basic analysis module ``dnds_analysis.py`` is provided which includes the
basic functions described below.

Initializing the analysis module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Having performed a scan using an instance of ``nptfit.NPTF``, the first thing to do is load the scan parameters in. This is done with

.. code:: python

    >>> nptf.load_scan()

where ``nptf`` is an instance of ``nptfit.NPTF``. 

.. NOTE::
   The analysis can be performed on an already existing scan. To do this an instance of ``nptfit.NPTF`` should be created with the same template and model configuration as used to perform the scan. Then instead of running the scan, simply load it and proceed with the analysis as below.

An instance of the analysis module can then be created as follows:

.. code:: python

    >>> an = dnds_analysis.Analysis(nptf, mask=None, pixarea=0.)

where ``mask`` specifies the ROI used for the analysis if this is different from that used for the scan, and ``pixarea`` is the area of a pixel in sr
if using non-``HEALPix`` maps.

Making triangle plots
~~~~~~~~~~~~~~~~~~~~~

Triangle/corner plots can be used to visualize multidimensional samples using a scatterplot matrix. 
A triangle plot with the default options using the `corner` package can be made using:

.. code:: python

    >>> an.make_triangle()

To use your own custom plotting options, use ``corner`` as follows

.. code:: python

    >>> corner.corner(an.nptf.samples, labels=an.nptf.params, range=[1 for i in range(an.nptf.n_params)])

with additional arguments as specified in http://corner.readthedocs.io/en/latest/.

Getting template intensities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Template intensities [counts/cm\ :sup:`2`/s/sr] can be calculated with

.. code:: python

    >>> an.return_intensity_arrays_poiss(comp)
    >>> an.return_intensity_arrays_non_poiss(comp)

for Poissonian and non-Poissonian templates with the key ``comp`` respectively. This returns a list of intensities corresponding to the posterior parameters
for the given template.

The non-Poissonian templates (NPT) intensity is calculated by integrating up :math:`\int_{S_{min}}^{S_{max}} dS~S~dN/dS`. This is approximated as a sum between :math:`S_{min}` and :math:`S_{max}`. The options associated with the non-Poissonian template intensity are:

+--------------+--------------+--------------+
| Argument     | Default      | Purpose      |
+==============+==============+==============+
| ``comp``     | \-           | The NPT key  |
+--------------+--------------+--------------+
| ``smin``     | 0.01         | Minimum      |
|              |              | counts       |
|              |              | to sum       |
|              |              | up from      |
+--------------+--------------+--------------+
| ``smax``     | 10000        | Maximum      |
|              |              | counts       |
|              |              | to sum       |
|              |              | up to        |
+--------------+--------------+--------------+
| ``nsteps``   | 10000        | Number of    |
|              |              | bins in s    |
|              |              | while summing|
|              |              | up           | 
+--------------+--------------+--------------+


Getting source count distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The posterior arrays for the source count distributions :math:`dN/dF` [counts\ :sup:`-1` cm\ :sup:`2` s deg\ :sup:`-2`] associated with a given template ``comp`` at a given ``flux`` (in counts/cm\ :sup:`2`/s) can be obtained using

.. code:: python

    >>> an.return_dndf_arrays(comp,flux)

The source count distribution can be plotted with

.. code:: python

    >>> an.plot_source_count_median(comp, smin, smax, nsteps, spow, **kwargs)
    >>> an.plot_source_count_band(comp, smin, smax, nsteps, spow, qs, **kwargs)

The options being the same as for obtaining the NPT intensity above. Additionally, spow is the power :math:`n` in :math:`F^ndN/dF` to return while plotting, and qs is an array of quantiles for which to return the dN/dF band.

Plotting intensity fractions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intensity fractions (fraction of template intensity to total intensity) for Poissonian and non-Poissonian templates respectively can be plotting using

.. code:: python

    >>> an.plot_intensity_fraction_poiss(comp, bins, **kwargs)
    >>> an.plot_intensity_fraction_non_poiss(comp, bins, **kwargs)

where ``comp`` is the template key, ``bins`` is the number of bins between 0 and 100 and ``**kwargs`` specify plotting options.


Accessing posteriors
~~~~~~~~~~~~~~~~~~~~

While the posteriors can be accessed with ``nptf.samples`` (or ``an.nptf.samples``) as above, the following functions provide a useful interfact to access individual parameters:

.. code:: python

    >>> an.return_poiss_parameter_posteriors(comp)
    >>> an.return_non_poiss_parameter_posteriors(comp)

where ``comp`` is the (non-)Poissonian template key.

For Poissonian models, this returns a list of posterior normalizaion parameters for that model. For non-Poissonian models, this returns three arrays:

.. code:: python

    >>> A_non_poiss_post, n_non_poiss_post, Sb_non_poiss_post = an.return_non_poiss_parameter_posteriors(comp)

where 

- ``A_non_poiss_post`` is an array of non-Poissonian normalization parameter posteriors
- ``n_non_poiss_post`` is a 2-D array, each sub-array containing posteriors for a given slope parameter, starting from the highest to the lowest
- ``Sb_non_poiss_post`` is a 2-D array, each sub-array containing posteriors for a given break parameter, starting from the highest to the lowest

Getting Bayesian log-evidences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Bayesian log-evidence and associated error can be accessed as follows:

.. code:: python

    >>> l_be, l_be_err = an.get_log_evidence()
