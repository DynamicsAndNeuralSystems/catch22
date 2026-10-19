% make_matlab_reference  Compute reference values for all 22 catch22 features
% with the MATLAB implementations in Matlab/, for every time series in
% tests/reference/series/. Writes tests/reference/matlab_reference.csv, which
% tests/check_matlab_reference.py compares against the C implementation.
%
% Run from MATLAB (no toolboxes required):
%   run('tests/reference/make_matlab_reference.m')

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..', 'Matlab')));

features = {'DN_HistogramMode_5', 'DN_HistogramMode_10', ...
    'DN_OutlierInclude_p_001_mdrmd', 'DN_OutlierInclude_n_001_mdrmd', ...
    'CO_f1ecac', 'CO_FirstMin_ac', 'SP_Summaries_welch_rect_area_5_1', ...
    'SP_Summaries_welch_rect_centroid', 'FC_LocalSimple_mean3_stderr', ...
    'FC_LocalSimple_mean1_tauresrat', 'MD_hrv_classic_pnn40', ...
    'SB_BinaryStats_mean_longstretch1', 'SB_BinaryStats_diff_longstretch0', ...
    'SB_MotifThree_quantile_hh', 'CO_HistogramAMI_even_2_5', 'CO_trev_1_num', ...
    'IN_AutoMutualInfoStats_40_gaussian_fmmi', 'SB_TransitionMatrix_3ac_sumdiagcov', ...
    'PD_PeriodicityWang_th0_01', 'CO_Embed2_Dist_tau_d_expfit_meandiff', ...
    'SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1', 'SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1'};

files = dir(fullfile(here, 'series', '*.txt'));
fid = fopen(fullfile(here, 'matlab_reference.csv'), 'w');
fprintf(fid, 'series,feature,value\n');
for f = 1:numel(files)
    y = load(fullfile(files(f).folder, files(f).name));
    y = (y - mean(y)) / std(y); % features are computed on the z-scored series
    [~, name] = fileparts(files(f).name);
    for k = 1:numel(features)
        fprintf(fid, '%s,%s,%.17g\n', name, features{k}, feval(features{k}, y));
    end
end
fclose(fid);
fprintf('Wrote reference values for %u series.\n', numel(files));
