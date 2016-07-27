function a = do_lpc(x,fp,do_gain_norm,lptype)
%*****************************************************************
% Usage a = do_lpc(x,fp,do_gain_norm,lptype)
% Frequency Domain LP implementation
%   x:    input signal 
%   fp:    order of the predictor
% LPC is implemented using autocorrelations (derived from the inverse
% Fourier transform of the power spectrum of the signal or least squares.
%*****************************************************************
% Sriram Ganapathy
% Johns Hopkins University
% Jan 2012
%*****************************************************************

if nargin < 2
    error(' Not Enough Input Parameters');
end
if nargin < 3
    do_gain_norm =0;
    lptype =1;
elseif nargin < 4
    lptype = 1;
end
flen = length(x);
if lptype == 1         % Auto-correlation LPC
    % Compute signal transform
    X = fft(x,2^nextpow2(2*flen-1));
    % Derive the autocorrelations using power spectrum
    R = ifft((X.*conj(X)/flen));

    % Derive the LPC model
    [a,g2] = levinson(R,fp);

    % Poles per column
    a = a.';
elseif lptype == 2    % Least Squares LPC
    x= x(:);
    N = flen;
    % Construct the Sample Matrix
    R = toeplitz ( x(fp:N-1), x(fp:-1:1)');
    % Obtain the least squares solution
    a = geninv(R) * x(fp+1:N);
    % Error norm for gain normalization 
    g2 = norm ( x(fp+1:N) - R*a);
    a = [1 ; -a];    
else
    error ('Unknown LP type. Use 1 or 2');
end

% Apply gain normalization if required (useful for channel distortions)
if ~(do_gain_norm)
   g = sqrt(g2);
   a = a/g;
end

