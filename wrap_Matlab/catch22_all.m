function [featureValues, featureNamesLong, featureNamesShort, isMasked] = catch22_all(data, doCatch24, maskNaNs)
% catch22_all   Calculate all catch22 (or catch24) features from an input time series
%
% Compile first with mexAll (once).
%
% ---INPUTS
% data: a univariate time series (vector)
% doCatch24: set to true to also include the mean and standard deviation (a total
%               of 24 features). Default: true
% maskNaNs: set a feature that cannot be estimated (e.g., for a constant or very
%               short series) to a fixed value at the appropriate end of its range,
%               so that real numbers are always returned. NaN/Inf in data is never
%               masked. Default: true
% ---OUTPUTS:
% featureValues: computed feature values for the input data (column vector)
% featureNamesLong, featureNamesShort: the feature names (cell arrays)
% isMasked: true for each masked feature value

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(doCatch24)
    doCatch24 = true;
end
if nargin < 3 || isempty(maskNaNs)
    maskNaNs = true;
end

%-------------------------------------------------------------------------------
% Compute all features with the compiled C implementation
%-------------------------------------------------------------------------------
[featureValues, isMasked] = catch22_mex(double(data(:)), doCatch24, maskNaNs);

if nargout > 1
    [featureNamesLong, featureNamesShort] = GetAllFeatureNames(doCatch24);
end

end
