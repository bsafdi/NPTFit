Initializing a scan
-------------------

The class ``nptfit.NPTF`` is used to instantiate, configure and run a
template fit. An instance of ``nptfit.NPTF`` can be created as follows:

.. code:: python

   >>> from NPTFit import nptfit
   >>> nptf = nptfit.NPTF(tag, work_dir, psf_dir)

Description of arguments (all optional):

+--------------+--------------+--------------+
| Argument     | Default      | Purpose      |
+==============+==============+==============+
| ``tag``      |``"Untagged"``| Label        |
|              |              | associated   |
|              |              | with the     |
|              |              | scan         |
+--------------+--------------+--------------+
| ``work_dir`` | ``None``     | Location     |
|              |              | where all    |
|              |              | the output   |
|              |              | from the     |
|              |              | scan is      |
|              |              | stored       |
+--------------+--------------+--------------+
| ``psf_dir``  | ``None``     | Where to     |
|              |              | save PSF     |
|              |              | correction   |
|              |              | files        |
+--------------+--------------+--------------+

.. NOTE::
   The current directory for ``work_dir`` and ``work_dir/psf_dir`` for 
   ``psf_dir`` are used if these are not provided.

