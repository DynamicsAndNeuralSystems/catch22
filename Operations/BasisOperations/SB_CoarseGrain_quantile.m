function yth = SB_CoarseGrain_quantile(y,numGroups)
% SB_CoarseGrain   Coarse-grains a continuous time series to a discrete alphabet.
%
%---INPUTS:
% howtocg = 'quantile' puts an equal number into each bin, the method of coarse-graining
%
% numGroups, either specifies the size of the alphabet for 'quantile' and 'updown'
%       or sets the timedelay for the embedding subroutines

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

N = length(y); % length of the input sequence

% ------------------------------------------------------------------------------
% Do the coarse graining
% ------------------------------------------------------------------------------

% quantile
th = quantile(y,linspace(0,1,numGroups+1)); % thresholds for dividing the time series values
th(1) = th(1)-1; % this ensures the first point is included
% turn the time series into a set of numbers from 1:numGroups
yth = zeros(N,1);
for i = 1:numGroups
    yth(y > th(i) & y <= th(i+1)) = i;
end

if any(yth == 0)
    error('All values in the sequence were not assigned to a group')
end

end
