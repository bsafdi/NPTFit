Initializing a scan
-------------------

The class ``nptfit.NPTF`` is used to instantiate, configure and run a
template fit. An instance of ``nptfit.NPTF`` can be created as follows:

.. code:: python

   >>> from NPTFit import nptfit
   >>> nptf = nptfit.NPTF(tag, work_dir)

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

.. NOTE::
   The current directory for ``work_dir`` is used if not provided.

