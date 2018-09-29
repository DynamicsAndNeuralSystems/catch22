function out = PD_PeriodicityWang_th0_01(y_in)
% PD_PeriodicityWang    Periodicity extraction measure of Wang et al.
%
% Implements an idea based on the periodicity extraction measure proposed in:
%
% "Structure-based Statistical Features and Multivariate Time Series Clustering"
% Wang, X. and Wirth, A. and Wang, L.
% Seventh IEEE International Conference on Data Mining, 351--360 (2007)
% DOI: 10.1109/ICDM.2007.103
%
% Detrends the time series using a single-knot cubic regression spline
% and then computes autocorrelations up to one third of the length of
% the time series. The frequency is the first peak in the autocorrelation
% function satisfying a set of conditions.
% 
% This code uses the SPLINEFIT v1.13.0.0 function by Jonas Lundgren from
% Matlab Central accessed on 25 Aug 2018.
%
%---INPUT:
% y, the input time series.
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

% ------------------------------------------------------------------------------
%% Foreplay
% ------------------------------------------------------------------------------
N = length(y_in); % length of the time series
th = 0.01; % [0, 0.01, 0.1, 0.2, 1/sqrt(N), 5/sqrt(N), 10/sqrt(N)]; % the thresholds with which to count a peak

% ------------------------------------------------------------------------------
%% 1: Detrend using a regression spline with 3 knots
% ------------------------------------------------------------------------------

ppSpline = splinefit(y_in');

y_spl = ppval(ppSpline,1:N);
% x = (ppSpline.breaks(1):ppSpline.breaks(2))-1;
% y_spl_1 = ppSpline.coefs(1,1)*x.^3 + ppSpline.coefs(1,2)*x.^2 + ppSpline.coefs(1,3)*x.^1 + ppSpline.coefs(1,4);
% x = (ppSpline.breaks(2):ppSpline.breaks(3))-ppSpline.breaks(2)+1;
% y_spl_2 = ppSpline.coefs(2,1)*x.^3 + ppSpline.coefs(2,2)*x.^2 + ppSpline.coefs(2,3)*x.^1 + ppSpline.coefs(2,4);
% y_spl = [y_spl_1, y_spl_2];

y = y_in - y_spl';

%% 2. Compute autocorrelations up to 1/3 the length of the time series.
acmax = ceil(N/3); % compute autocorrelations up to this lag
acf = zeros(acmax,1); % the autocorrelation function
for i = 1:acmax % i is the \tau, the AC lag
    acf(i) = mean(y(1:N-i).*y(i+1:N));
end

% ------------------------------------------------------------------------------
%% 3. Frequency is the first peak satisfying the following conditions:
% ------------------------------------------------------------------------------
% (a) a trough before it
% (b) difference between peak and trough is at least 0.01
% (c) peak corresponds to positive correlation

% (i) find peaks and troughs in ACF
diffac = diff(acf); % differenced time series
sgndiffac = sign(diffac); % sign of differenced time series
bath = diff(sgndiffac); % differenced, signed, differenced time series
troughs = find(bath == 2) + 1; % finds troughs
peaks = find(bath == -2) + 1; % finds peaks
npeaks = length(peaks);

out = 0;

for i = 1:npeaks % search through all peaks for one that meets the condition
    ipeak = peaks(i); % acf lag at which there is a peak
    thepeak = acf(ipeak); % acf at the peak
    ftrough = find(troughs < ipeak,1,'last');
    if isempty(ftrough); continue; end
    itrough = troughs(ftrough); % acf lag at which there is a trough (the first one preceeding the peak)
    thetrough = acf(itrough); % acf at the trough

    % (a) a trough before it: should be implicit in the ftrough bit above
    %     if troughs(1)>ipeak % the first trough is after it
    %         continue
    %     end

    % (b) difference between peak and trough is at least 0.01
    if thepeak - thetrough < th
        continue
    end

    % (c) peak corresponds to positive correlation
    if thepeak < 0
        continue
    end

    % we made it! Use this frequency!
    out = ipeak; break
end


end
