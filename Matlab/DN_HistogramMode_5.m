function out = DN_HistogramMode_5(y)

% no combination of single functions
coder.inline('never');

out = DN_HistogramMode(y, 5);