%-------------------------------------------------------------------------------
% Get the function names
featureNames = GetAllFeatureNames();
numFeatures = length(featureNames);

%-------------------------------------------------------------------------------
% Get the data
dataFileNames = {'../testData/test.txt', '../testData/test2.txt'};
numTestFiles = length(dataFileNames);

fprintf(1,'Testing %u compiled features on %u data files\n',numFeatures,numTestFiles);

%-------------------------------------------------------------------------------
% Test all functions on all the data files
for j = 1:numTestFiles

    dataFileName = dataFileNames{j};
    fprintf('\n%s\n\n', dataFileName)

    fileID = fopen(dataFileName,'r');
    data = fscanf(fileID,'%f');

    for featureInd = 1:numFeatures

        featureName = featureNames{featureInd};
        fh = str2func(['catch22_', featureName]);
        out = fh(data');

        fprintf("%s: %1.6f\n", featureName, out);
    end

    fclose(fileID);

end
