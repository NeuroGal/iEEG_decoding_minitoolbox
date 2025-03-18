# iEEG_decoding_minitoolbox
Matlab code for _single subject time-resolved decoding with statistics_ - useful for intracranial EEG (iEEG\ECoG) or spike data when you have one 'super subject'.

The decoding is performed using the MVPA_Light toolbox (citation below, function ‘decoding_wrapper.m’), but the stats are all implemented here to adapt to do within subject statistics (not group level) - the main function for this is ‘decoding_stats.m’ which calls many functions from my repository **time_resolved_stats** (https://github.com/NeuroGal/time_resolved_stats).

MVPA Light: Treder, MS Frontiers in Neuroscience (2020) https://doi.org/10.3389/fnins.2020.00289. Code here: https://github.com/treder/MVPA-Light (worked with a version downloaded on 04/02/2021).

This was written as part of the analysis for the following manuscript (so please cite if you use it): Vishne et al., Cell Reports 2023, 'Distinct Ventral Stream and Prefrontal Cortex Representational Dynamics during Sustained Conscious Visual Perception': https://doi.org/10.1016/j.celrep.2023.112752.

If you want to see how it was used in the paper visit https://github.com/NeuroGal/PersistentViewing_paper (folder: ‘3. Categories (decoding)’) - this also includes functions for plotting the output (see ‘Helper Functions’: nice_line_plot.m for time-point by time-point decoding and plot_results_mat.m for temporal generalization results).


Implemented statistics options (for more details see decoding_stats.m and **time_resolved_stats** repository):
- Analytic solution - relevant for accuracy results only - implemented by the function binomial_test.m in this repository. If you use this through decoding_stats.m false discovery rate is corrected using FDR.m (see more details below).
- Point-by-point permutations - implemented in the function perm_pvals.m in the time_resolved_stats repository. FDR correction when used through decoding_stats.m as in the analytic case.
- Max-statistic correction for multiple comparisons, using the function max_stat_correction.m in the time_resolved_stats repository. See Nichols & Holmes, 2002; https://doi.org/10.1002/hbm.1058 for more details.
- Cluster-based permutations correction for multiple comparisons, using the function run_cluster_perm.m in the time_resolved_stats repository, see Maris & Oostenveld, 2007; https://doi.org/10.1016/j.jneumeth.2007.03.024 for more details. Credit is due to Edden Gerber, an early version of this function is based on his code (see also https://edden-gerber.github.io/eeg_functions/).

For the analytic and point-by-point permutation cases FDR correction is done according to Benjamini and Hochberg (1995), using a function written by Edden Gerber (https://github.com/edden-gerber/time_series_analysis_and_statistics/blob/master/testing/FDR.m). I enclosed in the time_resolved_stats toolbox a version of this function which was on our lab server and definitely worked with this code.


Decoding is performed using MVPA_Light, called inside the function decoding_wrapper.m (which also calls the stats function):
- The function can receive multiclass data and either perform multiclass classification (if the requested classifier is multiclass_lda) or one vs one classification.
- The function runs by default temporal generalization (King and Dehaene, 2014; https://doi.org/10.1016/j.tics.2014.01.002), but also performs statistics for the diagonals separately in case you didn’t really want the full generalization matrix and you don’t want to do the correction for multiple comparisons for all of it.
- You can ask for multiple quantification metrics for the decoding (as implemented in MVPA_Light).
If this is all redundant for your usage simply call decoding_stats.m directly, but note you need to first generate the permutation data yourself.


Gal Vishne, June 2023

Twitter: @neuro_gal
