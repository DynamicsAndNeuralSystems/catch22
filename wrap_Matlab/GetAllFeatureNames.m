function [featureNamesLong,featureNamesShort] = GetAllFeatureNames(doCatch24)
% Retrieve all feature names from featureList.txt

if nargin < 1 || isempty(doCatch24)
    doCatch24 = true;
end
%-------------------------------------------------------------------------------
% Retrieve all (long) feature names from featureList.txt
fid = fopen('../featureList.txt','r');
formatSpec = '%s%s';
dataIn = textscan(fid,formatSpec,'CommentStyle','#','EndOfLine','\r\n','CollectOutput',true);
fclose(fid);
dataIn = dataIn{1}; % Collect one big matrix of cells
featureNamesLong = dataIn(:,1);
featureNamesShort = dataIn(:,2);

%-------------------------------------------------------------------------------
% Filter out the two extra features: 'DN_Mean' and 'DN_Spread_Std'
isMeanStd = strcmp(featureNamesLong,'DN_Mean') | strcmp(featureNamesLong,'DN_Spread_Std');
assert(sum(isMeanStd) == 2)
if ~doCatch24
    keepMe = ~isMeanStd;
    featureNamesLong = featureNamesLong(keepMe);
    featureNamesShort = featureNamesShort(keepMe);
    fprintf(1,'Using catch22.\n');
else
    fprintf(1,'Using catch24.\n');
end

end
