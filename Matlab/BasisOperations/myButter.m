function [Z_out, P_out] = myButter(n, W)
% Digital Butterworth filter, either 2 or 3 outputs
% Jan Simon, 2014, BSD licence
% See docs of BUTTER for input and output
% Fast hack with limited accuracy: Handle with care!
% Until n=15 the relative difference to Matlab's BUTTER is < 100*eps
V = tan(W * 1.5707963267948966);
Q = exp((1.5707963267948966i / n) * ((2 + n - 1):2:(3 * n - 1)));
nQ = length(Q);

Sg = V ^ nQ;
Sp = V * Q;
Sz = [];

% Bilinear transform:
P = (1 + Sp) ./ (1 - Sp);
Z = repmat(-1, size(P));
if isempty(Sz)
   G = real(Sg / prod(1 - Sp));
else
   G = real(Sg * prod(1 - Sz) / prod(1 - Sp));
   Z(1:length(Sz)) = (1 + Sz) ./ (1 - Sz);
end

% From Zeros, Poles and Gain to B (numerator) and A (denominator):
Z_out = G * real(poly(Z'));
P_out = real(poly(P));
