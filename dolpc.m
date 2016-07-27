function y = dolpc(x,order)
%y = dolpc(x,order)
%
% compute autoregressive model from spectral magnitude samples
%
% rows(x) = critical band
% col(x) = frame
%
% row(y) = lpc a_i coeffs, scaled by gain
% col(y) = frame
%

[nbands,nframes] = size(x);

if nargin < 2
  order = 8;
end

% Calculate autocorrelation 
r = real(ifft([x;x([(nbands-1):-1:2],:)]));
% First half only
r = r(1:nbands,:);

% Find LPC coeffs by durbin
[y,e] = levinson(r, order);

% Normalize each poly by gain
y = y'./repmat(e',(order+1),1);
