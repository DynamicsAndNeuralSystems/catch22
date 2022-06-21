function [featureValues, featureNames] = catch22_all(data, doCatch24)
% catch22_all   Calculate all catch22 (or catch24) features from an input time series
%
% ---INPUTS
% data: a univariate time series (row vector)
% catch24: set to true to also include mean and standard deviation (a total of 24
%               features)
% ---OUTPUTS:
% featureValues: computed feature values for the input data

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 2
    doCatch24 = true;
end

%-------------------------------------------------------------------------------
% Define the features:
%-------------------------------------------------------------------------------
featureNames = GetAllFeatureNames(doCatch24);

%-------------------------------------------------------------------------------
% Compute all features from their local compiled implementations
%-------------------------------------------------------------------------------
numFeatures = length(featureNames);
featureValues = zeros(numFeatures,1);

for featureInd = 1:numFeatures
    featureName = featureNames{featureInd};
    fh = str2func(['catch22_', featureName]);
    featureValues(featureInd) = fh(data');
end

end
