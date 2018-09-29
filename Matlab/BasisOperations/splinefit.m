function pp = splinefit(y)
%SPLINEFIT Fit a spline to noisy data.
%   PP = SPLINEFIT(X,Y,BREAKS) fits a piecewise cubic spline with breaks
%   (knots) BREAKS to the noisy data (X,Y). X is a vector and Y is a vector
%   or an ND array. If Y is an ND array, then X(j) and Y(:,...,:,j) are
%   matched. Use PPVAL to evaluate PP.
%
%   PP = SPLINEFIT(X,Y,P) where P is a positive integer interpolates the
%   breaks linearly from the sorted locations of X. P is the number of
%   spline pieces and P+1 is the number of breaks.
%
%   OPTIONAL INPUT
%   Argument places 4 to 8 are reserved for optional input.
%   These optional arguments can be given in any order:
%
%   PP = SPLINEFIT(...,'p') applies periodic boundary conditions to
%   the spline. The period length is MAX(BREAKS)-MIN(BREAKS).
%
%   PP = SPLINEFIT(...,'r') uses robust fitting to reduce the influence
%   from outlying data points. Three iterations of weighted least squares
%   are performed. Weights are computed from previous residuals.
%
%   PP = SPLINEFIT(...,BETA), where 0 < BETA < 1, sets the robust fitting
%   parameter BETA and activates robust fitting ('r' can be omitted).
%   Default is BETA = 1/2. BETA close to 0 gives all data equal weighting.
%   Increase BETA to reduce the influence from outlying data. BETA close
%   to 1 may cause instability or rank deficiency.
%
%   PP = SPLINEFIT(...,N) sets the spline order to N. Default is a cubic
%   spline with order N = 4. A spline with P pieces has P+N-1 degrees of
%   freedom. With periodic boundary conditions the degrees of freedom are
%   reduced to P.
%
%   PP = SPLINEFIT(...,CON) applies linear constraints to the spline.
%   CON is a structure with fields 'xc', 'yc' and 'cc':
%       'xc', x-locations (vector)
%       'yc', y-values (vector or ND array)
%       'cc', coefficients (matrix).
%
%   Constraints are linear combinations of derivatives of order 0 to N-2
%   according to
%
%     cc(1,j)*y(x) + cc(2,j)*y'(x) + ... = yc(:,...,:,j),  x = xc(j).
%
%   The maximum number of rows for 'cc' is N-1. If omitted or empty 'cc'
%   defaults to a single row of ones. Default for 'yc' is a zero array.
%
%   EXAMPLES
%
%       % Noisy data
%       x = linspace(0,2*pi,100);
%       y = sin(x) + 0.1*randn(size(x));
%       % Breaks
%       breaks = [0:5,2*pi];
%
%       % Fit a spline of order 5
%       pp = splinefit(x,y,breaks,5);
%
%       % Fit a spline of order 3 with periodic boundary conditions
%       pp = splinefit(x,y,breaks,3,'p');
%
%       % Constraints: y(0) = 0, y'(0) = 1 and y(3) + y"(3) = 0
%       xc = [0 0 3];
%       yc = [0 1 0];
%       cc = [1 0 1; 0 1 0; 0 0 1];
%       con = struct('xc',xc,'yc',yc,'cc',cc);
%
%       % Fit a cubic spline with 8 pieces and constraints
%       pp = splinefit(x,y,8,con);
%
%       % Fit a spline of order 6 with constraints and periodicity
%       pp = splinefit(x,y,breaks,con,6,'p');
%
%   See also SPLINE, PPVAL, PPDIFF, PPINT

%   Author: Jonas Lundgren <splinefit@gmail.com> 2010

%   2009-05-06  Original SPLINEFIT.
%   2010-06-23  New version of SPLINEFIT based on B-splines.
%   2010-09-01  Robust fitting scheme added.
%   2010-09-01  Support for data containing NaNs.
%   2011-07-01  Robust fitting parameter added.

%Copyright (c) 2017, Jonas Lundgren
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% inputs
x = 1:length(y);
breaks = [1, floor(length(y)/2), length(y)];
n = 4;
beta = 0;
dim = 1;

% Evaluate B-splines
[breaksB,coefsB,nB] = splinebase(breaks,n);
% % Make piecewise polynomial
% base = mkpp(breaksB,coefsB,nB);

coder.varsize('coefsB',[8 4], [0 0]);

% just don't use the struct...
if nargin==2, nB = 1; else nB = nB(:).'; end
dlk=numel(coefsB); l=length(breaksB)-1; dl=prod(nB)*l; k=fix(dlk/dl+100*eps);
% if (k<=0)||(dl*k~=dlk)
%    error(message('MATLAB:mkpp:PPNumberMismatchCoeffs',...
%        int2str(l),int2str(d),int2str(dlk)))
% end

breaksBase = reshape(breaksB,1,l+1);
coefsBase = coefsB; % reshape(coefsB,dl,k);
lBase = l;
kBase = k;
dBase = nB;

pieces = lBase;

% try this to get around unbounded-error in Coder.
coder.varsize('coefsBase',[8 4], [0 0]);

base = mkpp(breaksBase,coefsBase,dBase); % lBase,kBase,
A_base = ppval(base,x);

% Bin data
[junk,ibin] = histc(x,[-inf,breaks(2:end-1),inf]); %#ok

% Sparse system matrix
mx = numel(x);
ii = [ibin; ones(n-1,mx)];
ii = cumsum(ii,1);
jj = repmat(1:mx,n,1);

A = zeros(pieces+n-1,mx);
linIndices = sub2ind(size(A), ii, jj);
A(linIndices)=A_base;

% flipPoint = breaks(2)-1;
% coder.varsize('flipPoint',[1 1], [0 0]);
% A(1:end-1, 1:flipPoint) = A_base(:, 1:flipPoint);
% A(2:end, (flipPoint+1):end) = A_base(:, (flipPoint+1):end);

% Solve Min norm(u*A-y)
u = lsqsolve(A,y,beta);

% Compute polynomial coefficients
ii = [repmat(1:pieces,1,n); ones(n-1,n*pieces)];
ii = cumsum(ii,1);
jj = repmat(1:n*pieces,n,1);
% C = sparse(ii,jj,base.coefs,pieces+n-1,n*pieces);

C = zeros(pieces+n-1,n*pieces);
linIndices = sub2ind(size(C), ii, jj);
C(linIndices)=base.coefs;

coefs = u*C;
coefs = reshape(coefs,[],n);

% Make piecewise polynomial
pp = mkpp(breaks,coefs,dim);


%--------------------------------------------------------------------------
function [breaks0,coefsOut,n] = splinebase(breaks0,n)
%SPLINEBASE Generate B-spline base PP of order N for breaks BREAKS

assert (n == 4);

breaks0 = breaks0(:);     % Breaks
% breaks0 = breaks0';      % Initial breaks
h0 = diff(breaks0);       % Spacing
pieces0 = numel(h0);      % Number of pieces
deg = n - 1;            % Polynomial degree

% distances between knots to replicate points to the sides
hcopy = repmat(h0,ceil(deg/pieces0),1);

% to the left
hl = hcopy(end:-1:end-deg+1);
bl = breaks0(1) - cumsum(hl);
% and to the right
hr = hcopy(1:deg);
br = breaks0(end) + cumsum(hr);

% Add breaks
breaks = [bl(deg:-1:1); breaks0; br];
h = diff(breaks);
pieces = numel(h);

assert (pieces == 8);

% Initiate polynomial coefficients
coefs = zeros(n*pieces,n);
coefs(1:n:end,1) = 1;

% Expand h
ii = [1:pieces; ones(deg,pieces)];
ii = cumsum(ii,1);
ii = min(ii,pieces);
H = h(ii(:));

% Recursive generation of B-splines
for k = 2:n
    % Antiderivatives of splines
    for j = 1:k-1
        coefs(:,j) = coefs(:,j).*H/(k-j);
    end
    Q = sum(coefs,2);
    Q = reshape(Q,n,pieces);
    Q = cumsum(Q,1);
    c0 = [zeros(1,pieces); Q(1:deg,:)];
    coefs(:,k) = c0(:);
    % Normalize antiderivatives by max value
    fmax = repmat(Q(n,:),n,1);
    fmax = fmax(:);
    for j = 1:k
        coefs(:,j) = coefs(:,j)./fmax;
    end
    % Diff of adjacent antiderivatives
    coefs(1:end-deg,1:k) = coefs(1:end-deg,1:k) - coefs(n:end,1:k);
    coefs(1:n:end,k) = 0;
end

% Scale coefficients
scale = ones(size(H));
for k = 1:n-1
    scale = scale./H;
    coefs(:,n-k) = scale.*coefs(:,n-k);
end

% Reduce number of pieces
pieces = pieces - 2*deg;

% Sort coefficients by interval number
ii = [n*(1:pieces); deg*ones(deg,pieces)];
ii = cumsum(ii,1);
coefsOut = coefs(ii(:),:);

coder.varsize('coefsOut',[8 4], [0 0]);




%--------------------------------------------------------------------------
function u = lsqsolve(A,y,beta)
%LSQSOLVE Solve Min norm(u*A-y)


% Solution
u = y/A;

% % Robust fitting
% if beta > 0
%     [m,n] = size(y);
%     alpha = 0.5*beta/(1-beta)/m;
%     for k = 1:3
%         % Residual
%         r = u*A - y;
%         rr = r.*conj(r);
%         rrmean = sum(rr,2)/n;
%         rrmean(~rrmean) = 1;
%         rrhat = (alpha./rrmean)'*rr;
%         % Weights
%         w = exp(-rrhat);
%         spw = spdiags(w',0,n,n);
%         % Solve weighted problem
%         u = (y*spw)/(A*spw);
%     end
% end

