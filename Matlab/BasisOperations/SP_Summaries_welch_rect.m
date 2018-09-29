function out = SP_Summaries_welch_rect(y)
% SP_Summaries  Statistics of the power spectrum of a time series
%
% The estimation can be done using a periodogram, using the periodogram code in
% Matlab's Signal Processing Toolbox, or a fast fourier transform, implemented
% using Matlab's fft code.
%
%---INPUTS:
% y, the input time series
%
% psdMeth, the method of obtaining the spectrum from the signal:
%               (i) 'periodogram': periodogram
%               (ii) 'fft': fast fourier transform
%               (iii) 'welch': Welch's method
%
% windowType='rect', the window to use:
%               (i) 'boxcar'
%               (ii) 'rect'
%               (iii) 'bartlett'
%               (iv) 'hann'
%               (v) 'hamming'
%               (vi) 'none'
%
% nf, the number of frequency components to include, if
%           empty (default), it's approx length(y)
%
% dologabs, if 1, takes log amplitude of the signal before
%           transforming to the frequency domain.
%
% doPower, analyzes the power spectrum rather than amplitudes of a Fourier
%          transform
%
%---OUTPUTS:
% Statistics summarizing various properties of the spectrum,
% including its maximum, minimum, spread, correlation, centroid, area in certain
% (normalized) frequency bands, moments of the spectrum, Shannon spectral
% entropy, a spectral flatness measure, power-law fits, and the number of
% crossings of the spectrum at various amplitude thresholds.

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

% % ------------------------------------------------------------------------------
% %% Check that a Curve-Fitting Toolbox license is available:
% % ------------------------------------------------------------------------------
% BF_CheckToolbox('curve_fitting_toolbox')

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
% if size(y,2) > size(y,1)
%     y = y'; % Time series must be a column vector
% end
% if nargin < 2 || isempty(psdMeth)
%     psdMeth = 'fft'; % fft by default
% end
% if nargin < 3 || isempty(windowType)
%     windowType = 'hamming'; % Hamming window by default
% end
% if nargin < 4
%     nf = [];
% end
% if nargin < 5 || isempty(dologabs)
%     dologabs = false;
% end
% 
% if dologabs % a boolean
%     % Analyze the spectrum of logarithmic absolute deviations
%     y = log(abs(y));
% end

doPlot = false; % plot outputs
Ny = length(y); % time-series length

% default output
out.centroid = NaN;
out.area_5_1 = NaN;

%-------------------------------------------------------------------------------
% Set window (for periodogram and welch):
%-------------------------------------------------------------------------------

window = ones(Ny,1); % rectwin(Ny);

% ------------------------------------------------------------------------------
% Compute the Fourier Transform
% ------------------------------------------------------------------------------

% Welch power spectral density estimate:
Fs = 1; % sampling frequency
N = 2^nextpow2(Ny);
% [S, f] = pwelch(y,window,[],N,Fs);
% substitute with own Welch spectrum
[S_, f_] = welchy(y, N, Fs, window);

w_ = 2*pi*f_'; % angular frequency
S_ = S_/(2*pi); % adjust so that area remains normalized in angular frequency space

if ~any(isfinite(S_)) % no finite values in the power spectrum
    % This time series must be really weird -- return NaN (unsuitable operation)...
    % warning('NaN in power spectrum? A weird time series.');
    return
end

% Ensure both w and S are row vectors:
if size(S_,1) > size(S_,2)
    S = S_';
else
    S = S_;
end
if size(w_,1) > size(w_,2)
    w = w_';
else
    w = w_;
end

% if doPlot
%     figure('color','w')
%     plot(w,S,'.-k'); % plot the spectrum
%     % Area under S should sum to 1 if a power spectral density estimate:
%     title(sprintf('Area under psd curve = %.1f (= %.1f)',sum(S*(w(2)-w(1))),var(y)));
% end

N = length(S);
dw = w(2) - w(1); % spacing increment in w

% ------------------------------------------------------------------------------
% Shape of cumulative sum curve
% ------------------------------------------------------------------------------
csS = cumsum(S);

% this original was not well received by coder.
% f_frac_w_max = @(f) w(find(csS >= csS(end)*f,1,'first'));

% At what frequency is csS a fraction 0.5 of its maximum?
csSThres = csS(end)*0.5;
for i = 1:length(csS)
    if csS(i) > csSThres
        out.centroid = w(i);
        break
    end
end

% % wInd = find((csS >= csS(end)*0.5),1, 'first');
% % out = w(wInd); % f_frac_w_max(0.5);
% out = NaN;

% 5 bands
split = buffer_(S,floor(N/5));
if size(split,2) > 5, split = split(:,1:5); end
out.area_5_1 = sum(split(:,1))*dw;

end
