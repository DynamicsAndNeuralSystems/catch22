function out = MD_hrv_classic_pnn40(y)
% MD_hrv_classic    Classic heart rate variability (HRV) statistics.
%
% Typically assumes an NN/RR time series in units of seconds.
%
%---INPUTS:
% y, the input time series.
%
%---OUTPUTS:
%  pNN40
%  cf. "The pNNx files: re-examining a widely used heart rate variability
%           measure", J.E. Mietus et al., Heart 88(4) 378 (2002)
%
% Code is heavily derived from that provided by Max A. Little:
% http://www.maxlittle.net/

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

% Standard defaults
diffy = diff(y);
N = length(y); % time-series length

% ------------------------------------------------------------------------------
% Calculate pNNx percentage
% ------------------------------------------------------------------------------
% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
% strange to do this for a z-scored time series...
% pnntime = 20;

Dy = abs(diffy) * 1000;

% Anonymous function to do the PNNx calcualtion:
PNNxfn = @(x) sum(Dy > x)/(N-1);

% % exceed  = sum(Dy > pnntime);
% out.pnn5  = PNNxfn(5); %sum(Dy > 5)/(N-1); % proportion of difference magnitudes greater than 0.005*sigma
% out.pnn10 = PNNxfn(10); %sum(Dy > 10)/(N-1);
% out.pnn20 = PNNxfn(20); %sum(Dy > 20)/(N-1);
% out.pnn30 = PNNxfn(30); %sum(Dy > 30)/(N-1);
% pnn40
out = PNNxfn(40); %sum(Dy > 40)/(N-1);


end
