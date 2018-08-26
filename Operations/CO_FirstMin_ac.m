function out = CO_FirstMin_ac(y)
% CO_FirstMin  Time of first minimum in a given correlation function
%
%---INPUTS:
% y, the input time series
%
% Note that selecting 'ac' is unusual operation: standard operations are the
% first zero-crossing of the autocorrelation (as in CO_FirstZero), or the first
% minimum of the mutual information function ('mi').

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

N = length(y); % Time-series length

% ------------------------------------------------------------------------------
% Define the autocorrelation function
% ------------------------------------------------------------------------------

% Autocorrelation implemented as CO_AutoCorr
corrfn = @(x) CO_AutoCorr(y,x,'Fourier');

% ------------------------------------------------------------------------------
% Search for a minimum
% ------------------------------------------------------------------------------
% (Incrementally through time lags until a minimum is found)

autoCorr = zeros(N-1,1); % preallocate autocorrelation vector
for i = 1:N-1
    % Calculate the auto-correlation at this lag:
    autoCorr(i) = corrfn(i);

    % Hit a NaN before got to a minimum -- there is no minimum
    if isnan(autoCorr(i))
%         warning('No minimum in %s [[time series too short to find it?]]',minWhat)
        out = NaN; return
    end

    % We're at a minimum:
    if i==2 && (autoCorr(2) > autoCorr(1))
        % already increases at lag of 2 from lag of 1: a minimum (since ac(0) is maximal)
        out = 1;
        return
    elseif (i > 2) && (autoCorr(i-2) > autoCorr(i-1)) && (autoCorr(i-1) < autoCorr(i)) % minimum at previous i
        out = i-1; % I found the first minimum!
        return
    end
end

% Still decreasing -- no minimum was found after searching all across the time series:
out = N;

end
