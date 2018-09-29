function out = SB_MotifThree_quantile_hh(y)
% SB_MotifThree     Motifs in a coarse-graining of a time series to a 3-letter alphabet
%
% (As SB_MotifTwo but with a 3-letter alphabet)
%
%---INPUTS:
% y, time series to analyze
% cgHow, the coarse-graining method to use. 'quantile': equiprobable alphabet by time-series value
%
%---OUTPUTS:
% Statistics on words of length 1, 2, 3, and 4.

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
%% Set defaults:
%-------------------------------------------------------------------------------

% quantile graining
yt = SB_CoarseGrain_quantile(y,3);

N = length(yt); % Length of the symbolized sequence derived from the time series

% So we have a vector yt with entries \in {1,2,3}

% ------------------------------------------------------------------------------
%% Words of length 1
% ------------------------------------------------------------------------------
r1 = cell(3,1); % stores ranges as vectors
% out1 = zeros(3,1); % stores probabilities as doubles

for i = 1:3
	r1{i} = find(yt == i);
% 	out1(i) = length(r1{i})/N;
end


% ------------------------------------------------------------------------------
%% Words of length 2
% ------------------------------------------------------------------------------
% Make sure ranges are valid for looking at the next one
for i = 1:3
	if (~isempty(r1{i})) && (r1{i}(end) == N)
        r1{i} = r1{i}(1:end-1);
    end
end

r2 = cell(3,3);
out2 = zeros(3,3);
for i = 1:3
	for j = 1:3
		r2{i,j} = r1{i}(yt(r1{i}+1) == j);
		out2(i,j) = length(r2{i,j})/(N-1);
	end
end

% hh
out = f_entropy(out2); % entropy of this result

%-------------------------------------------------------------------------------
function h = f_entropy(x)
    % entropy of a set of counts, log(0)=0
    h = -sum(x(x > 0).*log(x(x > 0)));
end

end
