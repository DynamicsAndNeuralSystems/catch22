function out = DN_OutlierInclude_001_mdrmd(y,thresholdHow)
% DN_OutlierInclude     How statistics depend on distributional outliers.
%
% Measures a range of different statistics about the time series as more and
% more outliers are included in the calculation according to a specified rule,
% of outliers being furthest from the mean, greatest positive, or negative
% deviations.
%
% The threshold for including time-series data points in the analysis increases
% from zero to the maximum deviation, in increments of 0.01*sigma (by default),
% where sigma is the standard deviation of the time series.
%
% At each threshold, the mean, standard error, proportion of time series points
% included, median, and standard deviation are calculated, and outputs from the
% algorithm measure how these statistical quantities change as more extreme
% points are included in the calculation.
%
%---INPUTS:
% y, the input time series (ideally z-scored)
%
% thresholdHow, the method of how to determine outliers:
%     (i) 'abs': outliers are furthest from the mean,
%     (ii) 'p': outliers are the greatest positive deviations from the mean, or
%     (iii) 'n': outliers are the greatest negative deviations from the mean.
%
% inc, the increment to move through (fraction of std if input time series is
%       z-scored)
%
% Most of the outputs measure either exponential, i.e., f(x) = Aexp(Bx)+C, or
% linear, i.e., f(x) = Ax + B, fits to the sequence of statistics obtained in
% this way.
%
% [future: could compare differences in outputs obtained with 'p', 'n', and
%               'abs' -- could give an idea as to asymmetries/nonstationarities??]

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
%% Check Inputs
% ------------------------------------------------------------------------------
% If time series is all the same value -- ridiculous! ++BF 21/3/2010
if all(y == y(1)) % the whole time series is just a single value
%     fprintf(1,'The time series is a constant!\n');
    out = NaN; return % this method is not suitable for such time series: return a NaN
end

N = length(y); % length of the time series

inc = 0.01;

% ------------------------------------------------------------------------------
%% Initialize thresholds
% ------------------------------------------------------------------------------
switch thresholdHow

    case 'p' % analyze only positive deviations
        thr = (0:inc:max(y));
        tot = sum(y >= 0);

    otherwise % case 'n' % analyze only negative deviations
        thr = (0:inc:max(-y));
        tot = sum(y <= 0);
        
end

if isempty(thr)
    error('I suspect that this is a highly peculiar time series?!!!')
end

% msDt4 = zeros(length(thr),1); % median of index relative to middle
%                             
% for i = 1:length(thr)
%     th = thr(i); % the threshold
% 
%     % Construct a time series consisting of inter-event intervals for parts
%     % of the time serie exceeding the threshold, th
% 
%     if strcmp(thresholdHow,'n')% look at only positive deviations
%         r = find(y <= -th);
%     elseif strcmp(thresholdHow,'p')% look at only negative deviations
%         r = find(y >= th);
%     end
%     
%     msDt4 = median(r)/(N/2)-1;
% 
% end

msDt = zeros(length(thr),6); % mean, std, proportion_of_time_series_included,
                             % median of index relative to middle, mean,
                             % error
for i = 1:length(thr)
    th = thr(i); % the threshold

    % Construct a time series consisting of inter-event intervals for parts
    % of the time serie exceeding the threshold, th

    if strcmp(thresholdHow,'n')% look at only positive deviations
        r = find(y <= -th);
    else % if strcmp(thresholdHow,'p')% look at only negative deviations
        r = find(y >= th);
    end
    
    Dt_exc = diff(r); % Delta t (interval) time series; exceeding threshold

    msDt(i,1) = mean(Dt_exc); % the mean value of this sequence
    msDt(i,3) = length(Dt_exc)/tot*100; % this is just really measuring the distribution
                                      % : the proportion of possible values
                                      % that are actually used in
                                      % calculation
    msDt(i,4) = median(r)/(N/2)-1;


end

% ------------------------------------------------------------------------------
%% Trim
% ------------------------------------------------------------------------------
% Trim off where the number of events is only one; hence the differenced
% series returns NaN
%fbi = find(isnan(msDt(:,1)),1,'first'); % first bad index
fbi = findFirstTrue(isnan(msDt(:,1)));
if fbi ~= 0
    msDt = msDt(1:fbi-1,:);
    thr = thr(1:fbi-1);
end

% Trim off where the statistic power is lacking: less than 2% of data
% included
trimthr = 2; % percent
% mj = find(msDt(:,3) > trimthr,1,'last');
mj = findLastTrue(msDt(:,3) > trimthr);
if mj ~= 0
    msDt = msDt(1:mj,:);
    thr = thr(1:mj);
end

% ------------------------------------------------------------------------------
%%% Generate output:
% ------------------------------------------------------------------------------

%% Stationarity assumption

% mean, median and std of the median and mean of range indices
% out.mdrm = mean(msDt(:,4));
% out.mdrmd = median(msDt(:,4));
% out.mdrstd = std(msDt(:,4));
% 
% out.mrm = mean(msDt(:,5));
% out.mrmd = median(msDt(:,5));
% out.mrstd = std(msDt(:,5));

% mdrmd
out = median(msDt(:,4)); % median(msDt4(:));


end

function out = findFirstTrue(y)

    out = 0;
    
    coder.varsize('out', [1 1], [0 0]);

    for i = 1:length(y)
        if y(i)
            out = i;
            return
        end
    end

end

function out = findLastTrue(y)

    out = 0;
    
    coder.varsize('out', [1 1], [0 0]);

    for i = length(y):-1:1
        if y(i)
            out = i;
            return
        end
    end

end
