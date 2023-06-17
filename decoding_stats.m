function [mask, p_values, thresh] = decoding_stats(cfg_stats, main_stat, perm_stat)
% Statistics for *single subject* results from classification - useful for
%   intracranial EEG (iEEG\ECoG) or spike data when you have one 'super subject'.
% Can be used for checking difference from chance, or difference between
%   two decoding outputs. (in the second case insert as stat & perm_stat
%   the difference, so you should be doing perms for each case for each
%   case separately).
%
% Includes correction for multiple comparisons using cluster permutations 
% (Maris & Oostenveld, 2007; https://doi.org/10.1016/j.jneumeth.2007.03.024),
% using max-statistic control (Nichols & Holmes, 2002; https://doi.org/10.1002/hbm.1058),
% or correction for False Discovery Rate according to Benjamini and Hochberg (1995)
% with code by Edden Gerber (https://github.com/edden-gerber/time_series_analysis_and_statistics/blob/master/testing/FDR.m)
%
% Please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.
% 
% Input:
%   cfg_stats - Settings structure, see fields below.
%   main_stat - The decoding result of your real data.
%                   Size: n_time1 x n_time2 (can be n_time x 1 if it's not a temporal generalization result)
%   perm_stat - Permutation results - generate them before this function!
%                   Size: n_time1 x n_time2 x n_perm 
%                       * if you are using it for a vector don't forget to inset perms as
%                         n_time1 x 1 x n_perm (or 1 x ntime2 x perm, according to your main_stat)
%
% Output:
%   mask     - n_time1 x n_time2 (if you asked for cluster permutations and multiple 
%                                 clusters then this is incorporating all of them together)
%   p_values - n_time1 x n_time2 or array of cluster p_values (ordered from
%                                 largest to smallest cluster)
%   thresh   - what threshold was used in the FDR correction (in this case
%              the mask is created by comparing the p-value array to this 
%              threshold but the p-value array stays the same).
% 
%
% cfg_stats fields:
%   acc          - Default: false.
%                  True\False - says if binomial (analytic) is relevant.
%                  (use false if it's diff of accuracies, it's only for comparing to chance.
%                  relevant for both regualr (2-class) classification and multiclass and when
%                  the input is part of a confusion matrix)
%                  also determines the firstlevel calculation if you chose
%                  cluster-based permutations.
%   stat_type    - Default: 'analytic' if possible, otherwise 'perm'.
%                  Four options implemented:
%                   'analytic' - uses binomial test (only if possible:
%                                accuracy\parts of confusion matrices)
%                   'perm' - point by point comparison to the permutation
%                 -> both corrected for multiple comparisons using FDR correction
%                   'cluster' - runs cluster permutations (needs more inputs, see below)
%                   'max' - max statistic control for multiple comparisons
%   p_thresh     - Default: 0.05
%                  * This is for the cluster\max threshold or fdr, not the cluster first level. 
%   n_sides      - Default: 1 (right side)
%                  Options: 1 or 2. Decides if we check only higher from 
%                   chance (right sided test) or double sided.
%   chance_level - if acc == true: Default: 0.5 (Note this is not true for
%                  multiclass or for difference of accuracies!!!). Adjust accordingly!
%                  Also relevant for AUC, but then there is no default and you must insert it.  
% Analytic inputs (only relevant if acc == True):
%   n            - Needed for the binomial test, you get it from MVPA_Light.
%                  Must be user input (no default) [how many stimuli you
%                  had do construct the accuracy matrix]
% Cluster inputs:
%   n_clusters        - Default: Inf (all clusters). Max # of clusters to return. 
%                        In the permutations only the largest cluster is taken, but we 
%                        can choose to compare more than just our largest cluster
%                        (it's just stricter, but valid, see the original paper for more details)
%   firstlevel_type   - Default: 'p' if possible (accuracy\confusion, things
%                        the analytic solution can apply to). Otherwise: 'stat'
%                       How should the firstlevel_thresh be treated?
%                       Options: 'p'\'stat'\'zscore'
%                           'p' - p_value (only relevant here for things a
%                              binomial test can apply to)
%                           'stat' - apply the threshold to the decoding results 
%                              directly, e.g. anything above a specific auc level
%                           'zscore' - zscore the values according to the
%                              permutations mean & std and apply the
%                              threshold according to that. This is only
%                              applied in the function run_cluster_perm
%                              (and there it also setting chance_level to 0)
%   firstlevel_thresh - If p-value is possible, default: 0.05.
%                       Otherwise: MUST have user input.
%                       Note for firstlevel_type == 'stat':
%                       * If chance_level is given insert this as DISTANCE
%                         FROM CHANCE. (e.g. if you want all auc>=0.6 to be
%                         included and you inserted chance_level = 0.5 use
%                         0.1 as the firstlevel threshold).
%                       firstlevel_type == 'stat'\'zscore':
%                       * If n_sides == 2 this could be 2-entry array
%                         [pos, neg] <-POSITIVE FIRST, neg really has to
%                         be negative we're taking it relative to chance_level.
%                         If only one given it's applied as [+thresh,-thresh]
%                         (in both cases it's around chance_level if that is given)
%
% Written by Gal Vishne, Deouell Lab 2019-2021
% Bug reports \ requests: gal.vishne@gmail.com

cfg_stats = cfg_stats_parser(cfg_stats); % update defaults etc.

switch cfg_stats.stat_type
    case 'analytic'
        p_values = binomial_test(main_stat, cfg_stats);
    case 'perm'
        p_values = perm_pvals(main_stat, perm_stat, cfg_stats);
    case 'max'
        [mask, p_values, thresh] = max_stat_correction(main_stat, perm_stat, cfg_stats);
    case 'cluster'
        inputs = {cfg_stats, main_stat, perm_stat};
        if cfg_stats.acc && strcmp(cfg_stats.firstlevel_type,'p')
            p_values = binomial_test(main_stat, cfg_stats);
            perm_p_values = binomial_test(perm_stat, cfg_stats);
            inputs = [inputs, {p_values, perm_p_values}]; 
        end
        [cluster_masks, cluster_pvalues] = run_cluster_perm(inputs{:});
        cluster_masks = cat(3, cluster_masks{:}); % you can also export this if needed
        mask = any(cluster_masks, 3); 
        p_values = cluster_pvalues;
        thresh = [];
    otherwise
        error("No stat_type '%s'. Use 'analytic', 'perm', 'max' or 'cluster'", cfg_stats.stat_type)
end

% FDR correction & masks for the perm & analytic cases
if ismember(cfg_stats.stat_type, {'analytic','perm'})
    [~, thresh] = FDR(p_values(:), cfg_stats.p_thresh); % uses all values, change here if you want something else (e.g. take only above the diagonal since matrix is symmetric)    
    mask = p_values<=thresh;
end
end

function cfg_stats = cfg_stats_parser(cfg_stats)
% Make sure all fields are there & correct. Add default values for missing things.
if ~isfield(cfg_stats,'acc')
    warning('Assuming this is NOT an accuracy array'); cfg_stats.acc = false;
end    
if ~isfield(cfg_stats,'stat_type') || ~any(ismember(cfg_stats.stat_type,{'analytic','perm','cluster','max'}))
    if cfg_stats.acc
        warning('Setting stat_type to analytic'); cfg_stats.stat_type = 'analytic';
    else
        warning('Setting stat_type to perm'); cfg_stats.stat_type = 'perm';
    end
end
if ~cfg_stats.acc && strcmp(cfg_stats.stat_type, 'analytic')
    warning('No analytic test, setting stat_type to perm'); cfg_stats.stat_type = 'perm';
end
if cfg_stats.acc
    if strcmp(cfg_stats.stat_type, 'analytic') && ~isfield(cfg_stats,'n')
        error('To use binomial test you must set cfg_stats.n');
    end
    if ~isfield(cfg_stats,'chance_level')
        warning('Setting chance_level to 0.5'); cfg_stats.chance_level = 0.5;
    end  
end
if ~isfield(cfg_stats,'p_thresh')
    warning('Setting p_thresh to 0.05'); cfg_stats.p_thresh = 0.05;
end
if ~isfield(cfg_stats,'n_sides')
    warning('Setting n_sides to 1 (larger than, right sided only)'); cfg_stats.n_sides = 1;
end
if strcmp(cfg_stats.stat_type,'cluster')
    if ~isfield(cfg_stats,'n_clusters')
        warning('Setting n_clusters to Inf (all ordered from largest to smallest)'); cfg_stats.n_clusters = Inf;
    end
    if ~isfield(cfg_stats,'firstlevel_type') || ~any(ismember(cfg_stats.firstlevel_type,{'stat','p','zscore'}))
        if cfg_stats.acc
            warning('Setting firstlevel_type to p-value'); cfg_stats.firstlevel_type = 'p';
        else
            warning('Setting firstlevel_type to stat'); cfg_stats.firstlevel_type = 'stat';
        end
    end
    if ~isfield(cfg_stats,'firstlevel_thresh')
        if cfg_stats.acc
            warning('Setting firstlevel_thresh to 0.05'); cfg_stats.firstlevel_thresh = 0.05;
        else
            error('You must set field: firstlevel_thresh');
        end
    end
end
end