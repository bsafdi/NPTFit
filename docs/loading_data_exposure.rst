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

   >>> nptf.add_template(template,key)

where the first argument is the template ``numpy`` array and the second argument 
is the template key, used to identify the template in later calls.

All templates should be corrected for exposure before adding.  In detail, they should be models of the counts rather than the flux. Templates that will be used by Poissonian models should also be smoothed with the detector point spread function (PSF) in the case when that is is non-trivial.  Templates that will be used to describe the spatial distribution of non-Poissonian models should be exposure corrected but not smoothed with the PSF. 


.. WARNING::
   The data, exposure and template maps must be of the same length. 

.. NOTE::
   While the data, exposure and template maps need not be in 
   `HEALPix <http://healpix.jpl.nasa.gov/>`_ 
   format, PSF correction and creation of masks is only supported for HEALPix 
   formatted arrays.

.. _``HEALPix``: http://healpix.jpl.nasa.gov/
