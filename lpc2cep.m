function features = lpc2cep(a,nceps)
% features = lpc2cep(lpcas,nceps)
%    Convert the LPC 'a' coefficients in each column of lpcas
%    into frames of cepstra.
%    nceps is number of cepstra to produce, defaults to size(lpcas,1)

[nin, ncol] = size(a);

order = nin - 1;

if nargin < 2
  nceps = order + 1;
end

c = zeros(nceps, ncol);

% Code copied from HSigP.c: LPC2Cepstrum

% First cep is log(Error) from Durbin
c(1,:) = -log(a(1,:));

% Renormalize lpc A coeffs
a = a ./ repmat(a(1,:), nin, 1);
  
for n = 2:nceps
  sum = 0;
  for m = 2:n
    sum = sum + (n - m) * a(m,:) .* c(n - m + 1, :);
  end
  c(n,:) = -(a(n,:) + sum / (n-1));
end

features = c;

