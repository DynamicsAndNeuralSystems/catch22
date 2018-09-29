function out = IN_AutoMutualInfoStats_40_gaussian_fmmi(y)
% IN_AutoMutualInfoStats  Statistics on automutual information function for a time series.
%
%---INPUTS:
% y, column vector of time series data
%
% maxTau, maximal time delay
%
% estMethod = 'gaussian'
%
%---OUTPUTS:
% out = fmmi, statistic on the AMIs and their pattern across
%       the range of specified time delays.

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% no combination of single functions
coder.inline('never');

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------

% maxTau: the maximum time delay to investigate
maxTau = 40;

% Don't go above N/2
maxTau = min(maxTau,ceil(N/2));

% ------------------------------------------------------------------------------
%% Get the AMI data:
% ------------------------------------------------------------------------------

% do this manually, Gaussian really is the easiest MI
autocorrs = autocorrelation(y, maxTau);
ami = -0.5*log(1 - autocorrs.^2);


% ------------------------------------------------------------------------------
% Output statistics:
% ------------------------------------------------------------------------------
lami = length(ami);

% First minimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end-1).*dami(2:end) < 0 & dami(1:end-1)<0);
if isempty(extremai)
    % fmmi
    out = lami; % actually represents lag, because indexes don't but diff delays by 1
else
    % fmmi
    out = min(extremai);
end

end

function out = autocorrelation(y, maxLag)

    N = length(y);

    if maxLag >= N
        maxLag = N-1;
    end
    
    out = zeros(1,maxLag);
    
    for lag = 1:maxLag
        corrs = corrcoef(y(1:end-lag), y(1+lag:end));
        out(lag) = corrs(2,1);
    end

end
