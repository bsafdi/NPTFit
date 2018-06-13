Adding masks
------------

Masks can be used to restrict an analysis within a given region of interest 
(ROI). Masks are boolean arrays where pixels labelled as true are masked and 
those labelled false are unmasked. Masks can be loaded into an instance of 
``nptfit.NPTF`` (see :doc:`initializing_scan`) as:

.. code:: python

   >>> nptf.load_mask(mask)

.. WARNING::
   The mask must be the same length as the data.

Creating masks
-------------

For the case where the data is a 
`HEALPix <http://healpix.jpl.nasa.gov/>`_
array, masks can be created by the ``create_mask`` module. For example:

.. code:: python

    >>> from NPTFit import create_mask as cm
    >>> mask = cm.make_mask_total(nside, mask_ring = True, outer = 2, ring_b = -45)

where `nside` is the HEALPix resolution parameter. The following mask options are available. When set to ``True``, the associated 
arguments determine the parameters of the mask. If multiple options are set to
``True`` then the total mask is the combination of each option.

+-----------------+--------------+---------------+
| Argument        | Default      | Purpose       |
+=================+==============+===============+
| ``band_mask``   | ``False``    | If ``True``,  |
|                 |              | masks within  |
|                 |              | \|\ *b*\ \|   |
|                 |              | <             |
|                 |              | ``band_mask_  |
|                 |              | range``       |
+-----------------+--------------+---------------+
| ``l_mask``      | ``False``    | If ``True``,  |
|                 |              | masks         |
|                 |              | longitude     |
|                 |              | outside       |
|                 |              |               |
|                 |              | ``l_deg_min`` |
|                 |              | < *l* <       |
|                 |              | ``l_deg_max`` |
|                 |              |               |
+-----------------+--------------+---------------+
| ``b_mask``      | ``False``    | If ``True``,  |
|                 |              | masks         |
|                 |              | latitude      |
|                 |              | outside       |
|                 |              |               |
|                 |              | ``b_deg_min`` |
|                 |              | < *b* <       |
|                 |              | ``b_deg_max`` |
|                 |              |               |
+-----------------+--------------+---------------+
| ``mask_ring``   | ``False``    | If ``True``,  |
|                 |              | masks         |
|                 |              | outside       |
|                 |              | ``inner`` <   |
|                 |              | *r* <         |
|                 |              | ``outer``,    |
|                 |              |               |
|                 |              | of a ring     |
|                 |              | centred at    |
|                 |              | (``ring_b``,  |
|                 |              | \ ``ring_l``  |
|                 |              | )             |
+-----------------+--------------+---------------+
| ``custom_mask`` | ``None``     | Optional      |
|                 |              | user-provided |
|                 |              | mask          |
+-----------------+--------------+---------------+

.. NOTE::
   By convention, ``True`` pixels are masked and ``False`` unmasked.

.. NOTE::
   Creation of masks by ``create_mask`` is only supported for HEALPix 
   formatted maps. For other inputs, masks otherwise made can be used as 
   usual with ``nptf.load_mask()``.

.. TIP::
   See the :doc:`Example2_Creating_Masks` for an exposition of the options described here.


.. _``HEALPix``: http://healpix.jpl.nasa.gov/
