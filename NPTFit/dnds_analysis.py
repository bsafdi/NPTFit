###############################################################################
# dnds_analysis.py
###############################################################################
#
# Analyze results of a non-Poissonian template fit. Code to produce:
#
# - Template intensities and confidence intervals
# - Source count distributions
# - Intensity fractions
# - Triangle plots and log-evidences
#
# NB: code default to assumption analysis is done on a HEALPix map, if this
# is not the case, must insert a pixarea at initialization.
#
###############################################################################

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import corner


class Analysis:
    """ Class to analyze results of an NPTF.

        :param nptf: an instance of nptfit.NPTF, where load_scan() has been performed
        :param mask: if analysis is to be performed in a different ROI to the run, insert
        the analysis mask here
        :param pixarea: if using a non-HEALPix map, insert the area of a pixel (in sr)
    """

    def __init__(self, nptf, mask=None, pixarea=0.):

        self.nptf = nptf
        # Default to HEALPix map
        if pixarea == 0.:
            self.pixarea = 4*np.pi/self.nptf.npix
        else:
            self.pixarea = pixarea
        if mask is None:
            self.mask_total = self.nptf.mask_total
        else:
            self.mask_total = mask
        self.mask_compress_data()

    def return_intensity_arrays_non_poiss(self, comp, smin=0.01, smax=10000,
                                          nsteps=10000, counts=False):
        """ Return intensity quantiles of a non-Poissonian template

            :param comp: key of non-Poissonian template
            :param smin: minimum count to "integrate" dnds from
            :param smax: maximum count to "integrate" dnds to
            :param nsteps: number of count bins in sum approximation of integral
            :param counts: whether to return counts (or intensities, by default)
        """

        template = self.nptf.templates_dict_nested[comp]['template']
        template_masked_compressed = self.mask_and_compress(
                                            template, self.mask_total)

        # If intensity, convert from counts to counts/cm^2/s/sr
        if counts:
            self.template_sum = np.sum(template_masked_compressed)
        else:
            self.template_sum = np.mean(template_masked_compressed /
                                        self.exp_masked_compressed /
                                        self.pixarea)

        self.sarray = 10**np.linspace(np.log10(smin), np.log10(smax), nsteps)
        self.ds = [self.sarray[i+1]-self.sarray[i]
                   for i in range(len(self.sarray)-1)]
        self.ds = np.array(self.ds + [self.ds[-1]])

        # Get NPT intensity arrays. These are calculated as
        # \int(dS*S*dN/dS). Note that the APS parameter is a
        # rescaling of the counts, which is why to get the
        # intensity this is multiplied by the total counts
        self.intensity_array_non_poiss = \
            list(map(lambda sample:
                     np.sum(self.template_sum *
                            self.dnds(comp, sample, self.sarray) *
                            self.sarray*self.ds), self.nptf.samples))

        return self.intensity_array_non_poiss

    def return_intensity_arrays_poiss(self, comp, counts=False):
        """ Return intensity arrays of a Poissonian template

            :param comp: key of Poissonian template
            :param counts: whether to return counts (or intensities, by default)
        """

        template = self.nptf.templates_dict_nested[comp]['template']
        template_masked_compressed = self.mask_and_compress(
                                            template, self.mask_total)

        # If intensity, convert from counts to counts/cm^2/s/sr
        if counts:
            self.template_sum = np.sum(template_masked_compressed)
        else:
            self.template_sum = np.mean(template_masked_compressed /
                                        self.exp_masked_compressed /
                                        self.pixarea)

        # Get PT intensities by scaling the compressed mask intensity
        # by the relevant normalizations from chains
        self.intensity_array_poiss = \
            list(map(lambda sample: self.template_sum *
                     self.return_poiss_samples(comp, sample),
                     self.nptf.samples))

        return self.intensity_array_poiss

    def return_non_poiss_samples(self, comp, sample):
        """ Return non-Poissonian samples corrected for log priors
        """

        # Load all NPT models (stored after PT models)
        self.model_decompression_non_poiss = \
            np.array(self.nptf.model_decompression_key[self.nptf.n_poiss:])

        model_where = \
            np.where(self.model_decompression_non_poiss[:, 0] == comp)[0] \
            + self.nptf.n_poiss

        is_log_prior = \
            np.array(self.nptf.model_decompression_key)[model_where][:, 1]

        is_log_prior = list(is_log_prior == 'True')

        samples_model_not_log = \
            self.log_to_normal(np.array(sample)[model_where], is_log_prior)

        return samples_model_not_log

    def return_poiss_samples(self, comp, sample):
        """ Return Poissonian samples corrected for log priors
        """

        # Load all PT models
        self.model_decompression_poiss = \
            np.array(self.nptf.model_decompression_key[:self.nptf.n_poiss])

        model_where = \
            np.where(self.model_decompression_poiss[:, 0] == comp)[0]

        is_log_prior = \
            np.array(self.model_decompression_poiss)[model_where][0][1]

        samples_model_not_log = self.log_to_normal(
            np.array(sample)[model_where], [is_log_prior == 'True'])[0]

        return samples_model_not_log

    def return_dndf_arrays(self, comp, flux):
        """ Calcualte and return array of dN/dF values for the template comp
            and the given array of flux values (in counts/cm^2/s)
        """

        template = self.nptf.templates_dict_nested[comp]['template']
        template_masked_compressed = \
            self.mask_and_compress(template, self.mask_total)

        self.template_sum = np.sum(template_masked_compressed)

        # Rescaling factor to convert dN/dS to [(counts/cm^2 /s)^-2 /deg^2]
        # Note that self.area_mask has units deg^2.
        rf = self.template_sum*self.exp_masked_mean/self.area_mask

        # Get counts from flux
        s = np.array([flux])*self.exp_masked_mean

        return rf*np.array([self.dnds(comp, sample, s)[0]
                           for sample in self.nptf.samples])

    def calculate_dndf_arrays(self, comp, smin=0.01, smax=1000, nsteps=1000,
                              qs=[0.16, 0.5, 0.84]):
        """ Calculate dnds for specified quantiles
        """

        template = self.nptf.templates_dict_nested[comp]['template']
        template_masked_compressed = \
            self.mask_and_compress(template, self.mask_total)

        self.template_sum = np.sum(template_masked_compressed)

        self.sarray = 10**np.linspace(np.log10(smin), np.log10(smax), nsteps)

        self.flux_array = self.sarray/self.exp_masked_mean

        self.data_array = np.array([self.dnds(comp, sample, self.sarray)
                                   for sample in self.nptf.samples])

        # Rescaling factor to convert dN/dS to [(ph /cm^2 /s)^-2 /deg^2]
        # Note that self.area_mask has units deg^2.
        rf = self.template_sum*self.exp_masked_mean/self.area_mask

        self.qArray = [corner.quantile(self.data_array[::, i], qs)
                       for i in range(len(self.sarray))]

        self.qmean = rf*np.array([np.mean(self.data_array[::, i])
                                  for i in range(len(self.sarray))])

        self.qlow = rf*np.array([q[0] for q in self.qArray])
        self.qmid = rf*np.array([q[1] for q in self.qArray])
        self.qhigh = rf*np.array([q[2] for q in self.qArray])

    def plot_source_count_band(self, comp, smin=0.01, smax=1000, nsteps=1000,
                               qs=[0.16, 0.5, 0.84], spow=0, *args, **kwargs):
        """ Calculate and plot median source count function

            :param comp: key of non-Poissonian template
            :param smin: minimum count to plot
            :param smax: maximum count to plot
            :param nsteps: binning in counts s
            :param qs: source count quartles to plot
            :param spow: plotting s**spow*dn/ds
            **kwargs: plotting options
        """

        self.calculate_dndf_arrays(comp, smin=smin, smax=smax,
                                   nsteps=nsteps, qs=qs)
        plt.fill_between(self.flux_array, self.flux_array**spow*self.qlow,
                         self.flux_array**spow*self.qhigh, *args, **kwargs)

    def plot_source_count_median(self, comp, smin=0.01, smax=1000, nsteps=1000,
                                 spow=0, qs=[0.16, 0.5, 0.84], *args, **kwargs):
        """ Calculate and plot median source count function
        """

        self.calculate_dndf_arrays(comp, smin=smin, smax=smax,
                                   nsteps=nsteps, qs=qs)
        plt.plot(self.flux_array, self.flux_array**spow*self.qmid,
                 *args, **kwargs)

    def plot_intensity_fraction_non_poiss(self, comp, smin=0.00001, smax=1000,
                                          nsteps=1000, qs=[0.16, 0.5, 0.84],
                                          bins=50, color='blue',
                                          ls_vert='dashed', *args, **kwargs):
        """ Plot flux fraction of a non-Poissonian template

            :param bins: flux fraction bins
            :param color_vert: colour of vertical quartile lines
            :param ls_vert: matplotlib linestyle of vertical quartile lines
            **kwargs: plotting options
        """

        flux_fraction_array_non_poiss = \
            np.array(self.return_intensity_arrays_non_poiss(comp, smin=smin,
                     smax=smax, nsteps=nsteps, counts=True))/self.total_counts

        frac_hist_comp, bin_edges_comp = \
            np.histogram(100*np.array(flux_fraction_array_non_poiss), bins=bins,
                         range=(0, 100))

        qs_comp = \
            corner.quantile(100*np.array(flux_fraction_array_non_poiss), qs)

        plt.plot(bin_edges_comp[:-1],
                 frac_hist_comp/float(sum(frac_hist_comp)),
                 color=color, *args, **kwargs)

        for q in qs_comp:
            plt.axvline(q, ls=ls_vert, color=color)
        self.qs_comp = qs_comp

    def plot_intensity_fraction_poiss(self, comp, qs=[0.16, 0.5, 0.84], bins=50,
                                      color='blue', ls_vert='dashed',
                                      *args, **kwargs):
        """ Plot flux fraction of non-Poissonian component
        """

        flux_fraction_array_poiss = \
            np.array(self.return_intensity_arrays_poiss(comp, counts=True))\
            / self.total_counts

        frac_hist_comp, bin_edges_comp = \
            np.histogram(100*np.array(flux_fraction_array_poiss), bins=bins,
                         range=(0, 100))

        qs_comp = corner.quantile(100*np.array(flux_fraction_array_poiss), qs)

        plt.plot(bin_edges_comp[:-1],
                 frac_hist_comp/float(sum(frac_hist_comp)), 
                 color=color, *args, **kwargs)

        for q in qs_comp:
            plt.axvline(q, ls=ls_vert, color=color)
        self.qs_comp = qs_comp

    def return_poiss_parameter_posteriors(self, comp):
        """ Return posterior samples corresponding to individual parameters.

            :return: sample normalization values
        """

        self.samples_reduced_ary = [self.return_poiss_samples(comp, sample)
                                    for sample in self.nptf.samples]
        return self.samples_reduced_ary

    def return_non_poiss_parameter_posteriors(self, comp):
        """ Return posterior samples corresponding to individual parameters.

            :return: sample non-Poissonian posterior values values, list with
            self.aps_ary: sample normalization values
            self.n_ary_ary: sampled slopes, each sub-array corresponding to
                samples for that slope (highest to lowest).
            self.sb_ary_ary: sampled breaks, each sub-array corresponding to
                samples for that break (highest to lowest).
        """

        self.samples_reduced_ary = [self.return_non_poiss_samples(comp, sample)
                                    for sample in self.nptf.samples]
        self.samples_reduced_param_ary = list(zip(*self.samples_reduced_ary))

        nbreak = int((len(self.samples_reduced_ary[0]) - 2)/2.)

        self.aps_ary = self.samples_reduced_param_ary[0]

        self.n_ary_ary = [[] for _ in range(nbreak+1)]
        self.sb_ary_ary = [[] for _ in range(nbreak)]

        for i in range(nbreak+1):
            self.n_ary_ary[i] = self.samples_reduced_param_ary[i+1]

        for i in range(nbreak):
            self.sb_ary_ary[i] = self.samples_reduced_param_ary[i+nbreak+2]

        return self.aps_ary, self.n_ary_ary, self.sb_ary_ary

    def dnds(self, comp, sample, s):
        """ dN/dS values for NPT comp associated with a chain sample
        """

        samples_reduced = self.return_non_poiss_samples(comp, sample)

        nbreak = int((len(samples_reduced) - 2)/2.)

        # Get APS (float) and slopes/breaks (arrays)
        a_ps, n_ary, sb_ary = samples_reduced[0], samples_reduced[1:nbreak+2], \
            samples_reduced[nbreak+2:]

        # If relative breaks, define each break as (except the last one)
        # the multiplicative factor times the previous break
        if self.nptf.non_poiss_models[comp]['dnds_model'] \
                == 'specify_relative_breaks':
            for i in reversed(range(len(sb_ary) - 1)):
                sb_ary[i] = sb_ary[i+1]*sb_ary[i]

        # Determine where the s values fall with respect to the breaks
        where_vecs = [[] for _ in range(nbreak+1)]
        where_vecs[0] = np.where(s >= sb_ary[0])[0]
        for i in range(1, nbreak):
            where_vecs[i] = np.where((s >= sb_ary[i]) & (s < sb_ary[i-1]))[0]
        where_vecs[-1] = np.where(s < sb_ary[-1])[0]

        # Calculate dnds values for a broken power law with arbitrary breaks
        dnds = np.zeros(len(s))
        dnds[where_vecs[0]] = a_ps*(s[where_vecs[0]]/sb_ary[0])**(-n_ary[0])
        dnds[where_vecs[1]] = a_ps*(s[where_vecs[1]]/sb_ary[0])**(-n_ary[1])

        for i in range(2, nbreak+1):
            dnds[where_vecs[i]] = \
                a_ps*np.prod([(sb_ary[j+1]/sb_ary[j])**(-n_ary[j+1])
                              for j in range(0, i-1)]) * \
                (s[where_vecs[i]]/sb_ary[i-1])**(-n_ary[i])

        return dnds

    @staticmethod
    def log_to_normal(array, is_log):
        """ Take array and account for the impact of log priors
        """

        array_normal = []
        for i in range(len(array)):
            if is_log[i]:
                array_normal.append(10**array[i])
            else:
                array_normal.append(array[i])
        return array_normal

    @staticmethod
    def mask_and_compress(the_map, mask):
        """ Return compressed version of a map
        """

        map_masked = ma.masked_array(data=the_map, mask=mask)
        return map_masked.compressed()

    def mask_compress_data(self):
        """ Adjust the data and exposure for the mask
        """

        self.data_masked = ma.masked_array(data=self.nptf.count_map,
                                           mask=self.mask_total)
        self.data_masked_compressed = self.data_masked.compressed()
        self.total_counts = float(np.sum(self.data_masked_compressed))

        self.area_mask = np.sum(1-self.mask_total)*self.pixarea *\
            (360/(2.*np.pi))**2

        self.exp_masked = ma.masked_array(data=self.nptf.exposure_map,
                                          mask=self.mask_total)
        self.exp_masked_compressed = self.exp_masked.compressed()
        self.exp_masked_mean = np.mean(self.exp_masked_compressed)

    def get_log_evidence(self):
        """ Global log-evidence and associated error
        """

        self.lge = self.nptf.s['nested sampling global log-evidence']
        self.lge_err = self.nptf.s['nested sampling global log-evidence error']
        return self.lge, self.lge_err

    def make_triangle(self):
        """ Make a triangle plot
        """

        corner.corner(self.nptf.samples, labels=self.nptf.params, smooth=1.5,
                      smooth1d=1, quantiles=[0.16, 0.5, 0.84], show_titles=True,
                      title_fmt='.2f', title_args={'fontsize': 14},
                      range=[1 for _ in range(self.nptf.n_params)],
                      plot_datapoints=False, verbose=False)
