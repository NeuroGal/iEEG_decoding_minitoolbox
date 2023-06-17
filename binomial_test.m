function p_values = binomial_test(accs_array, cfg_stats)
% Compute binomial test p-values for single subject decoding - relevant for accuracy only.
%   * Important to leave the accuracies as proportion (out of 1) not move
%     to percent (out of 100).
%
% Inputs:
%   accs_array - your decoding accuracies for a single subject.
%                Size: n_time1 x n_time2 (n_time1\n_time2 can be == 1)
%   cfg_stats  - Settings structure, uses the following inputs (no defaults):
%                > n_sides - 1/2
%                   1-sided test is only the positive tail! (larger than)
%                   2-sided test checks both sides.
%                > chance_level (likely 0.5, but not always so you must 
%                   set it)
%                > n - number of items in the original classification (the
%                   rational is that if there is really no effect it is  
%                   harder to get a result very different from chance if 
%                   you have more items).
%
% Output: point-by-point p-values (no correction for multiple comparisons)
%
% Written by Gal Vishne, Deouell Lab ~2021
% Bug reports \ requests: gal.vishne@gmail.com
% Please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.

n = cfg_stats.n; chance_level = cfg_stats.chance_level;
accs_vec = floor(n*accs_array); E = n*chance_level;
if cfg_stats.n_sides == 1
    p_values = 1-binocdf(accs_vec-1, n, chance_level); % -1 to include probability of equal to (=)
elseif cfg_stats.n_sides == 2
    accs_abs = E+abs(accs_vec-E); % abs around chance
    p_values = 1-binocdf(accs_abs-1, n, chance_level);
end
end