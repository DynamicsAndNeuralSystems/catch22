% mexAll    Compile catch22_mex (used by catch22_all) from the C implementation.

% Path where the C implementation lives (relative to this file)
here = fileparts(mfilename('fullpath'));
basePath = fullfile(here, '..', 'C');

% All C files except the command-line tool
CFiles = dir(fullfile(basePath, '*.c'));
CFiles = CFiles(~strcmp({CFiles.name}, 'main.c'));
includeFiles = fullfile(basePath, {CFiles.name});

fprintf('Compiling catch22_mex from %u C files...\n', numel(includeFiles));
mex(['-I', basePath], '-outdir', here, fullfile(here, 'catch22_mex.c'), includeFiles{:});
fprintf(1, 'All done!\n');
