function [featureNamesLong, featureNamesShort, legacyNames, hctsaCalls] = GetAllFeatureNames(doCatch24)
% GetAllFeatureNames   The catch22 (or catch24) features, in the order of catch22_all
%
% ---INPUTS
% doCatch24: set to true to include the mean and standard deviation. Default: true
% ---OUTPUTS:
% featureNamesLong: long names (following hctsa's feature names)
% featureNamesShort: short names
% legacyNames: long names in catch22 v0.x (the definitions have since changed to
%               match hctsa, so values differ)
% hctsaCalls: the hctsa call that reproduces each feature exactly, on x_z (the
%               z-scored series) or x (the raw series)

if nargin < 1 || isempty(doCatch24)
    doCatch24 = true;
end

[featureNamesLong, featureNamesShort, legacyNames, hctsaCalls] = catch22_mex('features', doCatch24);

end
