% mexAll    Mex compile all .c files to be runnable from Matlab.

% Path where the C implementation lives
basePath = '../C';
ipath = {['-I', basePath], '-I.'};

% List all C files to include in the mex-call
CFiles = dir(basePath);
CFiles = CFiles(cellfun(@(x) contains(x, '.c'), {CFiles.name}));
CFileNames = {CFiles.name};

% Add path to the c file names
includeFiles = cellfun(@(x) fullfile(basePath, x), CFileNames, 'UniformOutput', false);

% Get function names
featureNames = GetAllFeatureNames(true);

% mex all feature functions separately
numFeatures = length(featureNames);
for i = 1:numFeatures
    featureName = featureNames{i};

    fprintf('[%u/%u]: Compiling %s...\n', i,numFeatures,featureName);
    mex(ipath{:}, ['catch22_', featureName,'.c'], 'M_wrapper.c', includeFiles{:})
    fprintf('\n');
end

fprintf(1,'All done!\n');
