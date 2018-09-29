function out = CO_AutoCorr(y,tau,whatMethod)
% CO_AutoCorr   Compute the autocorrelation of an input time series
%
%---INPUTS:
% y, a scalar time series column vector.
%
% tau, the time-delay. If tau is a scalar, returns autocorrelation for y at that
%       lag. If tau is a vector, returns autocorrelations for y at that set of
%       lags. Can set tau empty, [], to return the full function for the
%       'Fourier' estimation method.
%
% whatMethod, the method of computing the autocorrelation: 'Fourier',
%             'TimeDomainStat', or 'TimeDomain'.
%
%---OUTPUT: the autocorrelation at the given time-lag.
%
%---NOTES:
% Specifying whatMethod = 'TimeDomain' can tolerate NaN values in the time
% series.
%
% Computing mean/std across the full time series makes a significant difference
% for short time series, but can produce values outside [-1,+1]. The
% filtering-based method used by Matlab's autocorr, is probably the best for
% short time series, and is implemented here by specifying: whatMethod =
% 'Fourier'.

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
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2
    tau = 1; % Use a lag of 1 by default
end

if nargin < 3 || isempty(whatMethod)
    whatMethod = 'Fourier';
end

% ------------------------------------------------------------------------------
% Evaluate the time-series autocorrelation
% ------------------------------------------------------------------------------

switch whatMethod
case 'Fourier'
    % ------------------------------------------------------------------------------
    % Estimation based on Matlab function autocorr, based on method of:
    % [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
    %     Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
    %     NJ: Prentice-Hall, 1994.

    nFFT = 2^(nextpow2(length(y))+1);
    F = fft(y-mean(y),nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf./acf(1); % Normalize
    acf = real(acf);

    acf = acf(1:length(y));

    if isempty(tau) % return the full function
        out = acf;
    else % return a specific set of values
        out = zeros(length(tau),1);
        for i = 1:length(tau)
            if (tau(i) > length(acf)-1) || (tau(i) < 0)
                out(i) = NaN;
            else
                out(i) = acf(tau(i)+1);
            end
        end
    end



end

end
