% get function names
fid = fopen('../featureList.txt','r');
i = 1;
tline = fgetl(fid);
featureNames{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    featureNames{i} = tline;
end
fclose(fid);

% get the data
dataFileNames = {'../testData/test.txt', '../testData/test2.txt'};
for j = 1:2
    
    dataFileName = dataFileNames{j};
    fprintf('\n%s\n', dataFileName)
    
    fileID = fopen(dataFileName,'r');
    
    data = fscanf(fileID,'%f');
    
    for featureInd = 1:length(featureNames)-1
    
        featureName = featureNames{featureInd};

        fh = str2func(['catch22_', featureName]);

        out = fh(data');
        
        fprintf("%s: %1.6f\n", featureName, out);
    end
    
    fclose(fileID);
    
end
