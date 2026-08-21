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

% Compile the shared files (M_wrapper.c plus everything in basePath) to
% object files once, rather than recompiling all of them from source for
% every one of the ~24 features below:
sourcesToCache = [{'M_wrapper.c'}, includeFiles];
objDir = tempname;
mkdir(objDir);
cleanupObj = onCleanup(@() rmdir(objDir,'s'));
fprintf('Compiling %u shared C files (once)...\n',numel(sourcesToCache));
for i = 1:numel(sourcesToCache)
    mex(ipath{:}, '-c', '-outdir', objDir, sourcesToCache{i});
end
objFiles = dir(objDir);
objFiles = fullfile(objDir, {objFiles(~[objFiles.isdir]).name});

% Get function names
featureNames = GetAllFeatureNames(true);

% mex all feature functions separately, linking each against the cached
% shared objects instead of recompiling them
numFeatures = length(featureNames);
for i = 1:numFeatures
    featureName = featureNames{i};

    fprintf('[%u/%u]: Compiling %s...\n', i,numFeatures,featureName);
    mex(ipath{:}, ['catch22_', featureName,'.c'], objFiles{:})
    fprintf('\n');
end

fprintf(1,'All done!\n');
