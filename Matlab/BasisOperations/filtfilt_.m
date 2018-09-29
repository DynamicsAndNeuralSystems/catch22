function y = filtfilt_(b,a,x)
%FILTFILT Zero-phase forward and reverse digital IIR filtering.
%   Y = FILTFILT(B, A, X) filters the data in vector X with the filter
%   described by vectors A and B to create the filtered data Y.  The filter
%   is described by the difference equation:
%
%     a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%   The length of the input X must be more than three times the filter
%   order, defined as max(length(B)-1,length(A)-1).
%
%   References:
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed.,
%         McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-
%         backward filtering, IEEE Transactions on Signal Processing,
%         pp. 988-992, April 1996, Volume 44, Issue 4

%   Copyright 1988-2014 The MathWorks, Inc.


%% Parse coefficient coefficients vectors and determine initial conditions

na = numel(a);

L = 1;
% Check coefficients
b = b(:);
a = a(:);
nb = numel(b);
nfilt = max(nb,na);   
nfact = max(1,3*(nfilt-1));  % length of edge transients

% Zero pad shorter coefficient vector as needed
if nb < nfilt
    b(nfilt,1)=0;
elseif na < nfilt
    a(nfilt,1)=0;
end

% Compute initial conditions to remove DC offset
if nfilt>1
    zi = ( eye(nfilt-1) - [-a(2:nfilt), [eye(nfilt-2); ...
                                      zeros(1,nfilt-2)]] ) \ ...
     ( b(2:nfilt) - b(1)*a(2:nfilt) );
else
    zi = zeros(0,1);
end

%% Filter the data
for ii=1:L
    
    y = [2*x(1)-x(nfact+1:-1:2); x; 2*x(end)-x(end-1:-1:end-nfact)];
    
    % filter, reverse data, filter again, and reverse data again
    y = filter(b(:,ii),a(:,ii),y,zi(:,ii)*y(1));
    y = y(end:-1:1);
    y = filter(b(:,ii),a(:,ii),y,zi(:,ii)*y(1));
    
    % retain reversed central section of y
    y = y(end-nfact:-1:nfact+1);
end
