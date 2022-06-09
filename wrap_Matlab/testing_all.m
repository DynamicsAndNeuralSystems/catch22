% get the data
dataFileNames = {'../testData/test.txt', '../testData/test2.txt'};

for j = 1:2

    dataFileName = dataFileNames{j};
    fprintf('\n%s\n', dataFileName)

    fileID = fopen(dataFileName,'r');
    data = fscanf(fileID,'%f');
    fclose(fileID);

    % run all features
    [vals, names] = catch22_all(data);

    for featureInd = 1:length(names)
        fprintf("%s: %1.6f\n", names{featureInd}, vals(featureInd));
    end
end
