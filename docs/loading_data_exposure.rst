Loading data and exposure
-------------------------

All maps are passed as 1-d ``numpy`` arrays. The data
and exposure maps can be loaded into an instance of ``nptfit.NPTF`` (see :doc:`initializing_scan`) as
follows:

.. code:: python

   >>> nptf.load_data(data, exposure)

Adding templates
----------------

Spatial templates can be added as follows:

.. code:: python

   >>> nptf.add_template(template, key, units)

where the first argument is the template ``numpy`` array and the second argument 
is the template key, used to identify the template in later calls.

The argument `units` specifies the template units (counts or flux) or type (to be used in either Poissonian or PS models). The following values are allowed:

-- ``'counts'``: template in counts/pixel, to be used in a Poissonian model. Exposure and PSF corrected.
-- ``'flux'``: template in counts/cm^2/s/pixel, to be used in a Poissonian model. Not exposure corrected.
-- ``'PS'``: template for the underlying PS distribution, to be used in a non-Poissonian model. This shouldn't account for exposure effects. 

For example, a template specifying an underlying isotropic PS distribution for a non-Poissonian model would be added with the keyword ``'PS'`` as a truly isotropic array, e.g. ``[1,1,1,...,1]``.

.. WARNING::
   The data, exposure and template maps must be of the same length. 

.. NOTE::
   While the data, exposure and template maps need not be in 
   `HEALPix <http://healpix.jpl.nasa.gov/>`_ 
   format, PSF correction and creation of masks is only supported for HEALPix 
   formatted arrays.

.. _``HEALPix``: http://healpix.jpl.nasa.gov/
