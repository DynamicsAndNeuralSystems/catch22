function featureNames = GetAllFeatureNames(doCatch24)
% Retrieve all feature names from featureList.txt

if nargin < 1 || isempty(doCatch24)
    doCatch24 = true;
end
%-------------------------------------------------------------------------------

fid = fopen('../featureList.txt','r');
i = 1;
tline = fgetl(fid);
featureNames{i} = tline;
while ischar(tline)
    i = i + 1;
    tline = fgetl(fid);
    if ischar(tline)
        featureNames{i} = tline;
    end
end
fclose(fid);

%-------------------------------------------------------------------------------
% Filter out the two extra features: 'DN_Mean' and 'DN_Spread_Std'
isMeanStd = strcmp(featureNames,'DN_Mean') | strcmp(featureNames,'DN_Spread_Std');
assert(sum(isMeanStd) == 2)
if ~doCatch24
    featureNames = featureNames(~isMeanStd);
    fprintf(1,'Using catch22\n');
else
    fprintf(1,'Using catch24\n');
end

end
