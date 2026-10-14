function make_hctsa_reference(hctsaDir)
% make_hctsa_reference  Reference values for the catch24 features from hctsa.
%
% Evaluates each feature's hctsa call, as listed by `run_features --list` in
% tests/reference/features.tsv, on every series in tests/reference/series/, and
% writes tests/reference/hctsa_reference.csv (with the hctsa commit in
% hctsa_reference_version.txt). tests/check_matlab_reference.py compares the C
% (run with --no-mask, since hctsa returns NaN where catch22 masks) against it.
%
% Run from MATLAB with the path to an hctsa checkout:
%   make_hctsa_reference('/path/to/hctsa')
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(hctsaDir, 'Operations')), genpath(fullfile(hctsaDir, 'PeripheryFunctions')));

fid = fopen(fullfile(here, 'features.tsv'));
T = textscan(fid, '%s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fid);
names = T{1}; calls = T{4};

files = dir(fullfile(here, 'series', '*.txt'));
fid = fopen(fullfile(here, 'hctsa_reference.csv'), 'w');
fprintf(fid, 'series,feature,value\n');
for f = 1:numel(files)
    x = load(fullfile(files(f).folder, files(f).name)); x = x(:);
    x_z = (x - mean(x)) / std(x); %#ok<NASGU> (used by the calls)
    [~, series] = fileparts(files(f).name);
    for k = 1:numel(names)
        fprintf(fid, '%s,%s,%.17g\n', series, names{k}, SUB_eval(calls{k}, x, x_z));
    end
end
fclose(fid);

[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', hctsaDir));
if status ~= 0, commit = 'unknown'; end
fid = fopen(fullfile(here, 'hctsa_reference_version.txt'), 'w');
fprintf(fid, 'hctsa %s\nMATLAB %s\n', strtrim(commit), version);
fclose(fid);
fprintf('Wrote hctsa reference values for %u series.\n', numel(files));
end

function v = SUB_eval(call, x, x_z) %#ok<INUSD> (x and x_z are used by eval)
% A call is MATLAB code, optionally followed by .field to select an output.
tok = regexp(call, '^(.*\))\.(\w+)$', 'tokens', 'once');
if isempty(tok)
    v = eval(call);
else
    o = eval(tok{1});
    v = o.(tok{2});
end
end
