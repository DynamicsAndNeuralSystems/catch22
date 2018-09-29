function out = SB_BinaryStats_diff_longstretch0(y)
% SB_BinaryStats    Statistics on a binary symbolization of the time series
%
% Binary symbolization of the time series is a symbolic string of 0s and 1s.
%
% Provides information about the coarse-grained behavior of the time series
%
%---INPUTS:
% y, the input time series
%
% binaryMethod, the symbolization rule:
%         'diff': by whether incremental differences of the time series are
%                      positive (1), or negative (0),
%
%---OUTPUTS:
% Include the Shannon entropy of the string, the longest stretches of 0s
% or 1s, the mean length of consecutive 0s or 1s, and the spread of consecutive
% strings of 0s or 1s.

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

%-------------------------------------------------------------------------------
% Set defaults:
%-------------------------------------------------------------------------------

binaryMethod = 'diff';

%-------------------------------------------------------------------------------
% Binarize the time series:
%-------------------------------------------------------------------------------
yBin = BF_Binarize(y,binaryMethod);

%-------------------------------------------------------------------------------
% Consecutive string of ones / zeros (normalized by length)
%-------------------------------------------------------------------------------
difffy = diff(find([1;yBin;1]));
stretch0 = difffy(difffy ~= 1) - 1;

%-------------------------------------------------------------------------------
% pstretches
%-------------------------------------------------------------------------------

if isempty(stretch0) % all 1s (almost impossible to actually occur)
    longstretch0 = 0;
else
    longstretch0 = max(stretch0); % longest consecutive stretch of zeros
end

out = longstretch0;

end
