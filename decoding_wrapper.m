function [perf_all, stat_all, cat_order, perms_all] = decoding_wrapper(segs_array, class_vec, cfg_decoding, cfg_stats)
% Single subject decoding with stats - useful for intracranial EEG
%   (iEEG\ECoG) or spike data when you have one 'super subject'.
% Decoding uses MVPA_Light (citation below) - but the stats are different!
% Uses the function 'decoding_stats.m' that I wrote.
% 
% MVPA Light:
%   Treder, MS Frontiers in Neuroscience (2020) https://doi.org/10.3389/fnins.2020.00289
%   Code here: https://github.com/treder/MVPA-Light (worked with a version
%       downloaded on 04/02/2021)
% 
% Also please cite: Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% The code was written as part of the analysis for this paper.
%
% The function can recieve multiclass input and either performs one-vs-one decoding 
%   for all pairs or does multiclass decoding (depends on the input you use in
%   cfg_decoding). It automatically outputs both regular decoding and temporal
%   generalization (timextime) and if you ask for statistics it does both.
%   
% The code can also do cross decoding for the same stimuli in two different
%   repetitions of the same or different durations. This was specific to my
%   design and you probably don't need it...
% 
% Input:
%   - segs_array   - basic format: 3d array of n_stim x n_elec x n_time
%                     [other optional formats (cross decoding two
%                     repetitions of the same stim):
%                     same duration: 4d array n_stim x n_elec x n_time x 2
%                     different durations: cell array of 1x2 of the basic
%                     segment format]
%   - class_vec    - vector n_stim x 1 of the relevant classes - make
%                     sure numbers are 1-n_classes
%   - cfg_decoding - settings structure to be used by MVPA_Light (see the 
%                     options there), if the classifier is multiclass_lda 
%                     the code switches to multi_mode, otherwise one-vs-one
%                     decoding
%                     *note you can ask for multiple output metrics (as
%                     you can in MVPA_Light)* (# denoted as n_metric)
%   - cfg_stats (optional) - statistics settings structure, to be used in
%                     decoding_stats.m (after some adjustments; refer to
%                     that function for more details about the fields).
%                     > If no cfg_stats input is inserted no
%                       statistical quantification is performed.
%                     > When multiple output metrics are requested you can 
%                       insert a number of fields as cell arrays and the 
%                       content of each cell will be used to run the
%                       statistics for the corresponding decoding metric.
%                       Relevant fields: 'stat_type', 'firstlevel_type', 
%                       'firstlevel_thresh'.
%                     > The field n_perm defines how many permutations are
%                       run (which are given as input to decoding_stats, so
%                       it is not called there). Default: 1,000. (not used
%                       if only binomial test statistics are requested
%                       ('analytic' stat_type).
%                     > The following fields are adjusted inside the function
%                       before calling decoding_stats.m (in accordance with
%                       the settings requested in cfg_decoding): 'acc' (is
%                       accuracy data or not), 'n' (number of stimuli used
%                       in classification, relevant for analytic
%                       statistics), 'chance_level'.
%
% Output: 
%   - perf_all  - numeric array: n_time_1 x n_time_2 x n_comp x n_metric
%                  n_comp = # of comparisons
%                  if multi_mode this is a cell array (each metric separate
%                  cell since it can have a confusion matrix as output)
%   - stat_all  - cell: 2 (masks, p_values) x n_metric x n_comp x 2 (full time x time, diagonal)
%                  if multi_mode this won't have the n_comp dimension (and
%                  see below if you are asking for confusion matrix)
%   - cat_order - relevant for the regular mode, sort of legend: n_comp x 2
%                  (to know which two categories were compared in each case)
%   - perms_all - cell: n_metric x n_comp - I needed it for some things,
%                   not sure you will need it at any point
%
% Written by Gal Vishne, Deouell Lab 2021-2022
% Bug reports \ requests: gal.vishne@gmail.com

run_stats = exist('cfg_stats','var');
if (ndims(segs_array) == 4) % for cases of cross decoding, turn into cell format so the code after this can be just one version
    segs_array = squeeze(num2cell(segs_array, 1:3));
end
if iscell(segs_array) % the cross decoding case
    n_time = cellfun(@(x) size(x, 3), segs_array);
else % regular case
    n_time  = repmat(size(segs_array, 3),1,2);
end % in both cases we have 2 values now
if ~isfield(cfg_decoding, 'metric') % if no specific metric was requested use accuracy (the default in mvpa_light), needed for the stats
    cfg_decoding.metric = 'acc';
end
if isfield(cfg_decoding, 'classifier') && strcmp(cfg_decoding.classifier,'multiclass_lda')
    multi_mode = true;
else
    multi_mode = false;
end
if ~iscell(cfg_decoding.metric)
    n_metric = 1;
else
    n_metric = length(cfg_decoding.metric);
end
if run_stats
    if ~isfield(cfg_stats, 'n_perm')
            cfg_stats.n_perm = 1000;
            warning('Set n_perm to 1000 (not relevant for the analytic statistics option')
    end
    metrics = make_cell(cfg_decoding.metric, n_metric);
    stat_types = make_cell(cfg_stats.stat_type, n_metric);
    cluster_params = [];
    if any(strcmp(stat_types, 'cluster'))
        cluster_params.firstlevel_type = make_cell(cfg_stats.firstlevel_type, n_metric);
        cluster_params.firstlevel_thresh = make_cell(cfg_stats.firstlevel_thresh, n_metric);
    end
    acc_params = [];
end

if ~multi_mode
    unique_cat = unique(class_vec);
    cat_order  = nchoosek(unique_cat, 2);
    n_comp     = size(cat_order,1);
    perf_all   = nan(n_time(1), n_time(2), n_comp, n_metric);
    stat_all   = cell(2, n_metric, n_comp, 2);
    perms_all  = cell(n_metric, n_comp);
    remv_comp  = false(n_comp, 1); % if not enough stimuli to do a comparison
    for comp_num = 1:n_comp 
        cat1_ids = class_vec == cat_order(comp_num,1);
        cat2_ids = class_vec == cat_order(comp_num,2);
        if sum(cat1_ids)<2 || sum(cat2_ids)<2
            remv_comp(comp_num) = true;
            continue
        end
        relevant_trials = cat1_ids | cat2_ids;
        [X, clabel] = get_classif_inputs(segs_array, class_vec, relevant_trials);
        [perf, result_n] = inner_classif(X, clabel, cfg_decoding); % this is where we call MVPA_Light, see the inner function for what happens here
        if ~iscell(perf); perf = {perf}; end

        if run_stats
            % these are used if the metric is accuracy and we want an
            % analytic result (with binomial distribution)
            acc_params.n = result_n;
            acc_params.chance_level = 0.5; % NOTE: (I think) this depends on preprocessing that equates the class sizes
            update_inputs = {metrics, stat_types, acc_params, cluster_params}; % see the explanation in function update_cfg_stats below
            [mask, p_values, perms_all(:, comp_num)] = stats_wrapper(X, clabel, cfg_decoding, cfg_stats, update_inputs, perf);
            stat_all(1, :, comp_num, :) = mask; stat_all(2, :, comp_num, :) = p_values;
        else
            stat_all(:, :, comp_num, :) = repmat({nan(size(perf{1}))}, size(stat_all(:, :, comp_num, :)));
        end
        perf_all(:, :, comp_num, :) = cat(3, perf{:}); % not fit for confusion matrix (!)
        fprintf('Done with comparison %d/%d\n',comp_num,n_comp);
    end
    perf_all = perf_all(:, :, ~remv_comp, :); stat_all = stat_all(:,:,~remv_comp,:);
    perms_all = perms_all(:,~remv_comp); cat_order = cat_order(~remv_comp, :);
else % multi_mode
    [X, clabel] = get_classif_inputs(segs_array, class_vec, true(size(class_vec))); % last index - all trials are relevant, calling this to fix the class_vec numbering if needed
    [perf, result_n] = inner_classif(X, clabel, cfg_decoding);
    if ~iscell(perf); perf = {perf}; end
    perms_all  = cell(n_metric, 1);
    stat_all = repmat({nan(size(perf{1}))}, [2, n_metric, 2]); % fix size if not run_stats for downstream functions. won't fix for confusion
    
    if run_stats
        n_cat = length(unique(class_vec));
        acc_params.n = result_n; acc_params.chance_level = 1/n_cat; % NOTE: (I think) this depends on preprocessing that equates the class sizes
        update_inputs = {metrics, stat_types, acc_params, cluster_params}; % see the explanation in function update_cfg_stats below
        [mask, p_values, perms_all] = stats_wrapper(X, clabel, cfg_decoding, cfg_stats, update_inputs, perf);
        stat_all(1, :, :) = mask; stat_all(2, :, :) = p_values;
    end
    perf_all = perf; cat_order = [];
end

end

function [X, clabel] = get_classif_inputs(segs_array, class_vec, relevant_trials)
% get the segments and class labels for the current comparison
clabel = grp2idx(categorical(class_vec(relevant_trials))); % to turn it to 1,2 
if iscell(segs_array) % run_cross_decode
    X = cellfun(@(x) x(relevant_trials,:,:),segs_array,'UniformOutput',false);
else % run regular decode
    X = segs_array(relevant_trials, :, :);
end
end

function [perf, result_n] = inner_classif(X, clabel, cfg_decoding)
% here we call the mvpa_light functions - different inputs for the regular
% mode or the cross decoding mode (in both cases it does timextime)
if iscell(X) % run_cross_decode - IMPORTANT - we use each repetition once as train and once as test and average(!)
    [perf1, result] = mv_classify_timextime(cfg_decoding, X{1}, clabel, X{2}, clabel);
    perf2 = mv_classify_timextime(cfg_decoding, X{2}, clabel, X{1}, clabel);
    if iscell(perf1)
        perf = cell(size(perf1));
        for i = 1:length(perf1)
            perf{i} = (perf1{i} + perf2{i})/2;
        end
    else
        perf = (perf1 + perf2)/2;
    end
else % regular case
    [perf, result] = mv_classify_timextime(cfg_decoding, X, clabel);
end
result_n = result.n; % used in some cases for the stats later
end

function cfg_stats = update_cfg_stats(cfg_stats, idx, metrics, stat_types, acc_params, cluster_params)
% fixing up cfg_stats to be for one metric only (used inside the loop in stats_wrapper)
% idx - which metric are we doing now
cfg_stats.stat_type = stat_types{idx};
if strcmp(metrics{idx}, 'acc') || strcmp(metrics{idx}, 'confusion')
    cfg_stats.acc = true;
    cfg_stats.n = acc_params.n;
    cfg_stats.chance_level = acc_params.chance_level;
else
    cfg_stats.acc = false;
end
if strcmp(metrics{idx},'auc')
    cfg_stats.chance_level = acc_params.chance_level;
end
if ~(strcmp(metrics{idx},'acc') || strcmp(metrics{idx}, 'confusion')) && strcmp(cfg_stats.stat_type,'analytic')
    warning('The metric is not accuracy\confusion mat -> no analytic option to the test, moving to permutation stats (point by point with FDR correction)')
    cfg_stats.stat_type = 'perm';
end
if strcmp(cfg_stats.stat_type,'cluster')
    cfg_stats.firstlevel_thresh = cluster_params.firstlevel_thresh{idx};
    cfg_stats.firstlevel_type = cluster_params.firstlevel_type{idx};
end
end

function cell_version = make_cell(original_version, n_rep)
% for convenience
if iscell(original_version)
    cell_version = original_version;
else
    cell_version = repmat({original_version},n_rep,1);
end
end

function perm_stat = get_perms(X, clabel, cfg_decoding, n_perm, n_metric)
% shuffle the class labels to get permutations - takes a lot of time, but
% needed for stats in many cases. (calls to MVPA_Light)
perm_stat_interim = cell(n_metric, n_perm); % hack not to comit to dimensions in advance, merges to 1xn_metric cell array below
cfg_decoding.feedback = false; % so MVPA_Light won't spend so much time printing... 
% the CV status etc. You can change it, I print once every 25 repetitions below)

for p=1:n_perm
    perm_order = randperm(length(clabel)); % make sure your random seed works well! change here for complicated permutation schemes
    perm_stat_p = inner_classif(X, clabel(perm_order), cfg_decoding);
    if ~iscell(perm_stat_p); perm_stat_p = {perm_stat_p}; end
    perm_stat_interim(:,p) = perm_stat_p;
    if mod(p,25)==0; fprintf('Done with permutation %d/%d\n',p,n_perm); end
end
perm_stat = cell(1,n_metric);
rel_dim = cellfun(@ndims, perm_stat_interim(:,1)) + 1;
for met = 1:n_metric
    perm_stat{met} = cat(rel_dim(met), perm_stat_interim{met,:});
end
end

function [mask, p_values, perms] = stats_wrapper(X, clabel, cfg_decoding, cfg_stats, update_inputs, perf)
% The core happens in another function called decoding_stats. This mostly
%   generates permutation results if needed (by calling get_perms) and does
%   some messy stuff to account for the case of confusion matrix as output.
% The timextime matrix and diagonal are run separately so the diagonal
%   has a less strict correction for multiple comparisons if you didn't
%   really need the full matrix.
%
% update_inputs goes into update_cfg_stats: {metrics, stat_types, acc_params, cluster_params}
n_metric = length(update_inputs{1});
get_perms_cond = ... % only generate permutations if we need them - it takes a lot of time
    ~all(strcmp(update_inputs{2}, 'analytic')) || ~all(strcmp(update_inputs{1},'acc') | strcmp(update_inputs{1},'confusion'));
if get_perms_cond; perms = get_perms(X, clabel, cfg_decoding, cfg_stats.n_perm, n_metric); else; perms = cell(1,n_metric); end
mask = cell(n_metric, 2); p_values = cell(n_metric, 2); % 2nd dimension is (full matrix, diagonal only)
for met = 1:n_metric
    cfg_stats = update_cfg_stats(cfg_stats, met, update_inputs{:});
    stats_inputs = {cfg_stats, perf{met}}; stats_inputs{3} = perms{met}; % I don't remember why this is separate, but I assume there was some special case reason
	stats_inputs_diag = stats_inputs;
    if strcmp(update_inputs{1}{met},'confusion')
        conf_diags = cellfun(@diag, num2cell(permute(stats_inputs{2},[1 4 2 3]),[1 2]),'UniformOutput', false);
        stats_inputs_diag{2} = permute(cell2mat(conf_diags), [1 3 4 2]);
        conf_perm_diags = cellfun(@diag, num2cell(permute(stats_inputs{3},[1 4 2 3 5]),[1 2]),'UniformOutput', false);
        stats_inputs_diag{3} = permute(cell2mat(conf_perm_diags), [1 3 4 2 5]);
    else
        stats_inputs_diag{2} = diag(stats_inputs{2});
        perm_diags = cellfun(@diag, squeeze(num2cell(perms{met},[1 2])),'UniformOutput', false);
        stats_inputs_diag{3} = cat(3,perm_diags{:});
    end
    if strcmp(update_inputs{1}{met},'confusion') && strcmp(cfg_stats.stat_type,'cluster')
        % Special case to run each matrix separately (don't think of this as another dimention to connect clusters)
        % Notice currently we are not correcting for multiple comparisons across matrices!
        n_cat = size(perf{met},2); % n_categories
        mask_tmp_mat = nan(size(perf{met})); p_values_tmp_mat = cell(n_cat, n_cat);
        mask_tmp_diag = nan(size(stats_inputs_diag{2})); p_values_tmp_diag = cell(n_cat, n_cat);
        for cat1 = 1:n_cat
            for cat2 = 1:n_cat
                stats_inputs_conf = stats_inputs;
                stats_inputs_conf_diag = stats_inputs_diag;
                stats_inputs_conf{2} = squeeze(stats_inputs_conf{2}(:,cat1,cat2,:));
                stats_inputs_conf{3} = squeeze(stats_inputs_conf{3}(:,cat1,cat2,:,:));
                stats_inputs_conf_diag{2} = squeeze(stats_inputs_conf_diag{2}(:,cat1,cat2,:));
                stats_inputs_conf_diag{3} = permute(squeeze(stats_inputs_conf_diag{3}(:,cat1,cat2,:,:)),[1 3 2]);
                [mask_tmp_mat(:,cat1,cat2,:), p_values_tmp_mat{cat1,cat2}] = decoding_stats(stats_inputs_conf{:});
                [mask_tmp_diag(:,cat1,cat2,:), p_values_tmp_diag{cat1,cat2}] = decoding_stats(stats_inputs_conf_diag{:});
            end
        end
        mask{met,1} = mask_tmp_mat; mask{met,2} = mask_tmp_diag;
        p_values{met,1} = p_values_tmp_mat; p_values{met,2} = p_values_tmp_diag;
    else   
        [mask{met,1}, p_values{met,1}] = decoding_stats(stats_inputs{:}); % full matrix
        [mask{met,2}, p_values{met,2}] = decoding_stats(stats_inputs_diag{:}); % diagonal       
    end
end
end