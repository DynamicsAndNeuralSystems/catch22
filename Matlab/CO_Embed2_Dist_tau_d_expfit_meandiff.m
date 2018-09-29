function out = CO_Embed2_Dist_tau_d_expfit_meandiff(y)
% CO_Embed2_Dist    Analyzes distances in a 2-d embedding space of a time series.
%
% Returns statistics on the sequence of successive Euclidean distances between
% points in a two-dimensional time-delay embedding space with a given
% time-delay, tau.
%
% Outputs the the mean distance, the
% spread of distances, and statistics from an exponential fit to the
% distribution of distances.
%
%---INPUTS:
% y, a z-scored column vector representing the input time series.
% tau, the time delay.

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
%% Check inputs:
% ------------------------------------------------------------------------------

N = length(y); % time-series length

tau = CO_FirstZero(y,'ac');
if tau > N/10
    tau = floor(N/10);
end

% ------------------------------------------------------------------------------
% 2-dimensional time-delay embedding:
% ------------------------------------------------------------------------------

m = [y(1:end-tau), y(1+tau:end)];

% ------------------------------------------------------------------------------
% Calculate Euclidean distances between successive points in this space, d:
% ------------------------------------------------------------------------------

d = sqrt(diff(m(:,1)).^2 + diff(m(:,2)).^2); % Euclidean distance

% % ------------------------------------------------------------------------------
% % Empirical distance distribution often fits Exponential distribution quite well
% % Fit to all values (often some extreme outliers, but oh well)
% l = expfit(d);
l = mean(d); % maybe this is stupid, but I seriously don't see the difference

%% Sum of abs differences between exp fit and observed:

[N, binEdges] = histcounts_(d, [], true);

% from edges to centers
binCentres = mean([binEdges(1:end-1); binEdges(2:end)])';

% exponential fit in each bin
z = binCentres ./ l;
expf = exp(-z) ./ l;
expf(expf<0) = 0;

d_expfit_meandiff = mean(abs(N - expf)); % mean absolute error of fit

out = d_expfit_meandiff;

end
