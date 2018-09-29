function y_buff = buffer_(y, tau)

% replacement for y_buff = buffer(y,tau)

nCols = ceil(length(y)/tau);
yPadded = zeros(nCols*tau,1);
yPadded(1:length(y)) = y;
y_buff = reshape(yPadded, tau, nCols); 