function out = FC_LocalSimple_mean1_tauresrat(y)
% FC_LocalSimple    Simple local time-series forecasting.
%
% Simple predictors using the past trainLength values of the time series to
% predict its next value.
%
%---INPUTS:
% y, the input time series
%
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

% -------------------------------------------------------------------------
%% Set defaults
% -------------------------------------------------------------------------

trainLength = 1;

N = length(y); % Time-series length

% ------------------------------------------------------------------------------
% Do the local prediction
% ------------------------------------------------------------------------------
if strcmp(trainLength,'ac')
    lp = CO_FirstZero(y,'ac'); % make it tau
else
    lp = trainLength; % the length of the subsegment preceeding to use to predict the subsequent value
end
evalr = lp+1:N; % range over which to evaluate the forecast
if length(evalr)==0
%     warning('Time series too short for forecasting');
    out = NaN; return
end
res = zeros(length(evalr),1); % residuals
for i = 1:length(evalr)
    res(i) = mean(y(evalr(i)-lp:evalr(i)-1)) - y(evalr(i)); % prediction-value
end

% tauresrat
out = CO_FirstZero(res,'ac')/CO_FirstZero(y,'ac');

end
