function out = SC_FluctAnal_q2_taustep50_k1_logi_prop_r1(x,wtf,lag)
% SC_FluctAnal   Implements fluctuation analysis by a variety of methods.
%
% Much of our implementation is based on the well-explained discussion of
% scaling methods in:
% "Power spectrum and detrended fluctuation analysis: Application to daily
% temperatures" P. Talkner and R. O. Weber, Phys. Rev. E 62(1) 150 (2000)
%
% The main difference between algorithms for estimating scaling exponents amount
% to differences in how fluctuations, F, are quantified in time-series segments.
% Many alternatives are implemented in this function.
%
%---INPUTS:
% y, the input time series
%
% q, the parameter in the fluctuation function q = 2 (usual) gives RMS fluctuations.
%
% wtf, (what to fluctuate)
%       (i) 'endptdiff', calculates the differences in end points in each segment
%       (ii) 'range' calculates the range in each segment
%       (iii) 'std' takes the standard deviation in each segment
%
%           cf. "Evaluating scaled windowed variance methods for estimating the
%               Hurst coefficient of time series", M. J. Cannon et al. Physica A
%               241(3-4) 606 (1997)
%
%       (iv) 'iqr' takes the interquartile range in each segment
%       (v) 'dfa' removes a polynomial trend of order k in each segment,
%       (vi) 'rsrange' returns the range after removing a straight line fit
%
%           cf. "Analyzing exact fractal time series: evaluating dispersional
%           analysis and rescaled range methods",  D. C. Caccia et al., Physica
%           A 246(3-4) 609 (1997)
%
%       (vii) 'rsrangefit' fits a polynomial of order k and then returns the
%           range. The parameter q controls the order of fluctuations, for which
%           we mostly use the standard choice, q = 2, corresponding to root mean
%           square fluctuations.
%           An optional input parameter to this operation is a timelag for
%           computing the cumulative sum (or integrated profile), as suggested
%           by: "Using detrended fluctuation analysis for lagged correlation
%           analysis of nonstationary signals" J. Alvarez-Ramirez et al. Phys.
%           Rev. E 79(5) 057202 (2009)
%
% tauStep, increments in tau for linear range (i.e., if logInc = 0), or number of tau
%           steps in logarithmic range if login = 1
%           The spacing of time scales, tau, is commonly logarithmic through a range from
%           5 samples to a quarter of the length of the time series, as suggested in
%           "Statistical properties of DNA sequences", C.-K. Peng et al. Physica A
%           221(1-3) 180 (1995)
%
%           Max A. Little's fractal paper used L = 4 to L = N/2:
%           "Exploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder Detection"
%           M. A. Little et al. Biomed. Eng. Online 6(1) 23 (2007)
%
% k, polynomial order of detrending for 'dfa', 'rsrangefit'
%
% lag, optional time-lag, as in Alvarez-Ramirez (see (vii) above)
%
% logInc, whether to use logarithmic increments in tau (it should be logarithmic).
%
%---OUTPUTS: include statistics of fitting a linear function to a plot of log(F) as
% a function of log(tau), and for fitting two straight lines to the same data,
% choosing the split point at tau = tau_{split} as that which minimizes the
% combined fitting errors.
%
% This function can also be applied to the absolute deviations of the time
% series from its mean, and also for just the sign of deviations from the mean
% (i.e., converting the time series into a series of +1, when the time series is
% above its mean, and -1 when the time series is below its mean).
%
% All results are obtained with both linearly, and logarithmically-spaced time
% scales tau.

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

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(wtf)
    wtf = 'rsrange'; % re-scaled range analysis by default
end

if nargin < 3 || isempty(lag)
    lag = '';
end

q = 2;
tauStep = 50;
k = 1;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

N = length(x); % length of the time series
doPlot = 0; % plot relevant outputs to figure

% 1) Compute integrated sequence

if isempty(lag)
    % didn't specify a lag, do a normal cumsum:
    y = cumsum(x);
else
    % specified a lag, do a decimation:
    y = cumsum(x(1:lag:end));
end

%-------------------------------------------------------------------------------
% Perform scaling over a range of tau, up to a fifth the time-series length
%-------------------------------------------------------------------------------
% Peng (1995) suggests 5:N/4 for DFA
% Caccia suggested from 10 to (N-1)/2...
%-------------------------------------------------------------------------------
taur = unique(round(exp(linspace(log(5),log(floor(N/2)),tauStep))));
ntau = length(taur); % analyze the time series across this many timescales

if ntau < 8 % fewer than 8 points
%     fprintf(1,'This time series (N = %u) is too short to analyze using this fluctuation analysis\n',N);
    out = NaN;
    return
end

% 2) Compute the fluctuation function as follows
F = zeros(1,ntau);
% F is the fluctuation function
% each entry correponds to a given scale tau, and contains
% the fluctuation function at that scale

for i = 1:ntau
    % buffer the time series at the scale tau
    tau = taur(i); % the scale on which to compute fluctuations

    % y_buff = buffer(y,tau); % buff not supported by Coder
    nCols = ceil(length(y)/tau);
    yPadded = zeros(nCols*tau,1);
    yPadded(1:length(y)) = y;
    y_buff = reshape(yPadded, tau, nCols); 
    
    if size(y_buff,2) > floor(N/tau) % zero-padded, remove trailing set of points...
        y_buff = y_buff(:,1:end-1);
    end

    % analyzed length of time series (with trailing end-points removed)
    nn = size(y_buff,2)*tau;

    switch wtf
        case 'dfa'
            tt = (1:tau)'; % faux time range
            for j = 1:size(y_buff,2)
                % fit a polynomial of order k in each subsegment
                p = polyfit(tt,y_buff(:,j),k);
                % remove the trend, store back in y_buff
                y_buff(:,j) = y_buff(:,j) - polyval(p,tt);
            end
            % reshape to a column vector, y_dt (detrended)
            y_dt = reshape(y_buff,nn,1);
        case 'rsrangefit' % polynomial fit (order k) rather than endpoints fit: (~DFA)
            tt = (1:tau)'; % faux time range
            for j = 1:size(y_buff,2)
                % fit a polynomial of order k in each subsegment
                p = polyfit(tt,y_buff(:,j),k);
                % remove the trend, store back in y_buff
                y_buff(:,j) = y_buff(:,j) - polyval(p,tt);
            end
            y_dt = max(y_buff) - min(y_buff); % range(y_buff);
        otherwise
            error('Unknown fluctuation analysis method ''%s''',wtf);
    end

    % Compute fluctuation function:
    F(i) = (mean(y_dt.^q)).^(1/q);
end


%-------------------------------------------------------------------------------
% Smooth unevenly-distributed points in log space:
%-------------------------------------------------------------------------------
logtt = log(taur);
logFF = log(F);
ntt = ntau;

% ------------------------------------------------------------------------------
%% Try assuming two components (2 distinct scaling regimes)
% ------------------------------------------------------------------------------
% Move through, and fit a straight line to loglog before and after each point.
% Find point with the minimum sum of squared errors

% First spline interpolate to get an even sampling of the interval
% (currently, in the log scale, there are relatively more at large scales

% Deterine the errors
sserr = nan(ntt,1); % don't choose the end points
minPoints = 6;
for i = minPoints:ntt-minPoints
    r1 = 1:i;
    p1 = polyfit(logtt(r1),logFF(r1),1);
    r2 = i:ntt;
    p2 = polyfit(logtt(r2),logFF(r2),1);
    % Sum of errors from fitting lines to both segments:
    sserr(i) = norm(polyval(p1,logtt(r1))-logFF(r1)) + norm(polyval(p2,logtt(r2))-logFF(r2));
end

% breakPt is the point where it's best to fit a line before and another line after
breakPt = find(sserr == min(sserr),1,'first');
r1 = 1:breakPt;

prop_r1 = length(r1)/ntt;


% output
out = prop_r1;

end
