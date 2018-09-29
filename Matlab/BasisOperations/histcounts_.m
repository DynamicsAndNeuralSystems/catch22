function [N, binEdges] = histcounts_(d, nBinsIn, normalise)
    % own function for counting 

    if nargin<2 || isempty(nBinsIn)
        nBinsIn = [];
    end
    
    if nargin<3 || isempty(normalise)
        normalise = false;
    end

    
    dMinMax = [min(d), max(d)];
    dRange = dMinMax(2) - dMinMax(1);
    
    % fixed or automatic number of bins
    if isempty(nBinsIn)
        nBins = ceil(dRange/(3.5*std(d)/(numel(d)^(1/3))));
    else
        nBins = nBinsIn;
    end
    
    % compute edges
    binWidth = dRange/nBins;
    binEdges = dMinMax(1) + (0:nBins)*binWidth;
    
    % count occurences in bins
    N = zeros(nBins, 1);
    d_ = d - dMinMax(1);
    for i = 1:length(d_)
        binInd = floor(d_(i)/binWidth)+1;
        binInd = max(1, binInd);
        binInd = min(nBins, binInd);
        N(binInd) = N(binInd) + 1;
    end
    
    % probability normalisation
    if normalise
        N = N./numel(d);
    end