% testing   Check catch22_all against the C outputs in testData/basic (run after mexAll).

here = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', 'testData', 'basic');
dataFileNames = {'test', 'test2', 'testShort', 'testSinusoid'};
[featureNamesLong, featureNamesShort] = GetAllFeatureNames(true);
numFeatures = length(featureNamesLong);

fprintf(1, 'Testing %u features on %u data files.\n', numFeatures, length(dataFileNames));
numBad = 0;
for j = 1:length(dataFileNames)
    data = readmatrix(fullfile(dataDir, [dataFileNames{j}, '.txt']), 'FileType', 'text');
    featureValues = catch22_all(data, true);

    % Expected values: lines of "value, name, time" written by run_features --catch24
    expected = readcell(fullfile(dataDir, [dataFileNames{j}, '_output.txt']), ...
                        'FileType', 'text', 'Delimiter', ',');
    [isListed, ind] = ismember(featureNamesLong, strtrim(expected(:, 2)));
    assert(all(isListed), 'Not every feature is in %s_output.txt', dataFileNames{j});
    expectedValues = cell2mat(expected(ind, 1));

    isOK = abs(featureValues - expectedValues) <= 1e-9 * max(1, abs(expectedValues)) ...
            | (isnan(featureValues) & isnan(expectedValues));
    fprintf(1, '\n%s: %u/%u match the C\n', dataFileNames{j}, sum(isOK), numFeatures);
    for i = find(~isOK)'
        fprintf(1, '  %s (%s): %.14g (C: %.14g)\n', featureNamesShort{i}, ...
                featureNamesLong{i}, featureValues(i), expectedValues(i));
    end
    numBad = numBad + sum(~isOK);
end
assert(numBad == 0, '%u feature values differ from the C', numBad);
fprintf(1, '\nAll features match the C.\n');
