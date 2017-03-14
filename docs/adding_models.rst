Adding models
-------------

From the list of templates that ``nptfit.NPTF`` knows about (see :doc:`loading_data_exposure`), we can define
an arbitrary number of Poissonian (smooth/diffuse) and non-Poissonian (point source)
models.

Poissonian models
~~~~~~~~~~~~~~~~~

Poissonian models only have one parameter associated with them: the template normalisation.
Poissonian models corresponding to an added template can be loaded as:

.. code:: python

    nptf.add_poiss_model(template_name='iso', 
                         model_tag='$A_{iso}$', 
                         prior_range=[-2,2], 
                         log_prior=False, 
                         fixed=False, 
                         fixed_norm=1.0)

Arguments for Poissonian model:

+---------------------+---------------+---------------------------------------------------+
| Argument            | Default       | Purpose                                           |
+=====================+===============+===================================================+
| ``template_name``   |      \-       | Key corresponding to previously loaded template   |
+---------------------+---------------+---------------------------------------------------+
| ``model_tag``       |      \-       | ``LaTeX``-ready string, for plots                 |
+---------------------+---------------+---------------------------------------------------+
| ``prior_range``     | ``[min,max]`` |  Prior range to scan over                         |
+---------------------+---------------+---------------------------------------------------+
| ``log_prior``       | ``False``     | Whether to scan template in log-space             |
+---------------------+---------------+---------------------------------------------------+
| ``fixed``           | ``False``     | Whether the model is fixed                        |
+---------------------+---------------+---------------------------------------------------+
| ``fixed_norm``      | 1.0           | Template normalization if the model is fixed      |
+---------------------+---------------+---------------------------------------------------+

.. NOTE::
   When ``log_prior=True``, the associated values of ``prior_range`` and ``fixed_norm`` must also be the logs of the values being used.

.. WARNING::
   When ``log_prior=True``, the log used for the parameters is base 10.

Non-Poissonian models
~~~~~~~~~~~~~~~~~~~~~

The flux distribution of non-Poissonian (point source) models is modeled as a multiply broken 
power law with a specified number of breaks :math:`l`, the best-fit parameters of which can 
then be inferred. This source count distribution, which gives the differential number of sources per
unit of flux, takes the form

.. math:: 
    \frac{dN}{dF} = A \left\{ \begin{array}{lc} \left( \frac{F}{F_{b,1}} \right)^{-n_1}, & F \geq F_{b,1} \\ \left(\frac{F}{F_{b,1}}\right)^{-n_2}, & F_{b,1} > F \geq F_{b,2} \\ \left( \frac{F_{b,2}}{F_{b,1}} \right)^{-n_2} \left(\frac{F}{F_{b,2}}\right)^{-n_3}, & F_{b,2} > F \geq F_{b,3} \\ \left( \frac{F_{b,2}}{F_{b,1}} \right)^{-n_2} \left( \frac{F_{b,3}}{F_{b,2}} \right)^{-n_3} \left(\frac{F}{F_{b,3}}\right)^{-n_4}, & F_{b,3} > F \geq F_{b,4} \\ \\
    \ldots & \ldots \\ \\
    \left[ \prod_{i=1}^{\ell-1} \left( \frac{F_{b,i+1}}{F_{b,i}} \right)^{-n_{i+1}} \right] \left( \frac{F}{F_{b,\ell}} \right)^{-n_{\ell+1}} & F_{b,\ell} > F \end{array} \right.

It is sometimes convenient to specify the breaks in terms of counts instead of flux.  However, if the exposure map is non-uniform over the ROI, then the notion of counts in pixel dependent.  While the NPTFit code properly accounts for the pixel-dependent exposure correction, we also allow the user the specify the breaks :math:`F_{b,i}` in terms of an effective number of counts :math:`S_{b,i} \equiv F_{b,i} \cdot \text{mean}_\text{ROI}(E_p)`, where :math:`\text{mean}_\text{ROI}(E_p)` is the mean of the exposure map :math:`E_p` over the ROI.

Non-Poissonian models corresponding to an added template can be loaded as:

.. code:: python

    nptf.add_non_poiss_model(template_name='iso',
                             model_tag=['$A_{ps}$','$n_1$','$n_2$','$S_b^{(1)}$'], 
                             prior_range=[[-6,6],[2.05,30],[-2,1.95],[0.05,30.0]], 
                             log_prior=[True,False,False,False],
                             dnds_model='specify_breaks',
                             fixed_params=None,units='counts')

Arguments for non-Poissonian model:

+-------------------+--------------------+--------------------------------------------------------+
| Argument          | Default            | Purpose                                                |
+===================+====================+========================================================+
| ``template_name`` | \-                 | Key corresponding to loaded template                   |
+-------------------+--------------------+--------------------------------------------------------+
| ``model_tag``     | \-                 | LaTeX-ready string of nb-broken power-law parameters   |
|                   |                    |                                                        |
|                   |                    | [A, n_1, … , n_{nb+1}, Sb_1, … , Sb_nb] with Sb_1/n_1  |
|                   |                    |                                                        |
|                   |                    | the highest break/slope                                |
+-------------------+--------------------+--------------------------------------------------------+
| ``prior_range``   | ``[]``             | Prior range to scan over, given as a list of [min,max] |
|                   |                    |                                                        |
|                   |                    | each model parameter                                   |
+-------------------+--------------------+--------------------------------------------------------+
| ``log_prior``     | ``False``          | Whether to scan each parameter in log-space            |
+-------------------+--------------------+--------------------------------------------------------+
| ``dnds_model``    | ``specify_breaks`` | Whether to use absolute or relative breaks             |
|                   |                    |                                                        |
|                   |                    | (see :math:`dN/dS` model specifications below)         |
+-------------------+--------------------+--------------------------------------------------------+
| ``fixed_params``  | ``None``           | Which parameters to keep fixed (see Fixed parameter    | 
|                   |                    |                                                        |
|                   |                    | specifications below)                                  |
+-------------------+--------------------+--------------------------------------------------------+
| ``units``         |   ``counts``       | Whether the breaks in `dN/dF` are specified in terms of|
|                   |                    |                                                        |
|                   |                    | :math:`F_{b}` or :math:`S_{b}`, which is defnied above.|       
|                   |                    |                                                        |
|                   |                    | (see units specifications below)                       |
+-------------------+--------------------+--------------------------------------------------------+

.. NOTE::
   The number of breaks in the non-Poissonian model is inferred from the length of the ``model_tag`` array.

.. WARNING::
   Non-Poissonian (or PS) models must use a template loaded with ``units='PS'``, while non-Poissonian models should use ``units='counts'`` or ``units='flux'``.

:math:`dN/dF` model specifications
****************************

where :math:`n_1` is the highest index and :math:`S_b^{(1)}` the highest
break.


The following options are allowed for ``dnds_model``:

- ``specify_breaks``: all breaks are specified in absolute counts, :math:`\left[ A, n_1, \ldots, n_{\ell+1}, S_b^{(1)}, \ldots, S_b^{(\ell)} \right]`
- ``specify_relative_breaks``: the highest break is specified in counts, with each subsequent lower break specified relative to the subsequent higher break.
  :math:`\left[ A, n_1, \ldots, n_{\ell+1}, S_b^{(1)}, \lambda^{(2)}, \ldots, \lambda^{(\ell - 1)}, \lambda^{(\ell)} \right]` where :math:`\lambda^{(i)} = S_b^{(i)}/S_b^{(i-1)}`.

Fixed parameter specifications
******************************

Fixed parameters should be passed as an array with syntax 

.. code:: python

    fixed_params = [[param_index_1,fixed_value_1],[param_index_2,fixed_value_2]]

where parameter indexing starts from 0.

Units specifications
****************************

The following options are allowed for ``units``:

- ``counts``: The code assumes the user has specified the break in counts (:math:`S_b`) and will infer the breaks in flux (:math:`F_b`) by dividing the mean exposure: :math:`F_{b,i} \equiv S_{b,i} / \text{mean}_\text{ROI}(E_p)`. 
- ``flux``: The code assumes the breaks are already specified in terms of flux

.. TIP::
   See :doc:`Example5_Running_nonPoissonian_Scans` and :doc:`Example7_Galactic_Center_nonPoissonian` for examples of using these options in an analysis.
