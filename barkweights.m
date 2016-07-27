function [wts,idx] = barkweights(flen,sr,type,par,dB)
% [wts,idx] = barkweights(flen,sr,type,par)
%   flen: length of frame (dct or fft length)
%   sr:   sampling rate
%   type: 'plp' or 'gauss' for now
%   par:  parameter passed to the function, for now only used by the
%         gaussian for control of the width (meaningful values 1 wide, 2 steeper)
%   dB:   how many dB's down to set the threshold (default 0 keeps all) 
%
%   wts:  the weights in a cell array 
%   idx:  the indices to which the weights correspond (for speed up) 
%
%   notes: Dan Ellis' code



if nargin < 3; type = 'gauss'; end
if nargin < 4; par  = 1;       end
if nargin < 5; dB   = 48;      end

% This is in my case the DCT length
nfreqs = flen;

% How many output bands? (copied from rasta/init.c)
nyqbar = hz2bark(sr/2);
nfilts = ceil(nyqbar)+1;
% bark per filt
step_barks = nyqbar/(nfilts - 1);
% step_barks = 0.9734; %PM PRIDAT PRO 16kHz
% sr = 2*8.0009e+03; %PM PRIDAT PRO 16kHz

% How many bands to discount at each edge (usually 1)
first_good = round(1/step_barks);

% Bark frequency of every bin in FFT
binbarks = hz2bark([0:(nfreqs-1)]*(sr/2)/(nfreqs-1));

% Weights to collapse FFT bins into aud channels
wts = cell(nfilts,1);
% Initialize idx with full range
idx = repmat([1,flen],nfilts,1);

switch lower(type)
    case {'plp','dan','hynek'}
        for I = 1:nfilts
            f_bark_mid = (I-1) * step_barks;
            % Linear slopes in log-space (i.e. dB) intersect to trapezoidal window
            % half a bark to the left and half to the right [+0.5;-0.5]
            wts{I} = 10.^(min(0,min([binbarks - f_bark_mid + 0.5;-2.5*(binbarks - f_bark_mid - 0.5)])));
        end
    case {'gaus','gauss','gaussian'}
        N = flen;
        a = par(1); 
        % Index vector
        %k = -(N-1)/2:(N-1)/2;
        % Equation 44a from [1]
        %        wts = exp((-1/2)*(a * k/(N/2)).^2)';
        for I = 1:nfilts
            f_bark_mid = (I-1) * step_barks;
            wts{I} = exp((-1/2)*(a * (binbarks - f_bark_mid)).^2);
        end
    otherwise
        error('Unknown filter type, try "plp" or "gaussian"');
end

% Readjust windows to -dB threshold for speed
if dB
    % convert dB to linear scale
    lin = 10^(-dB/20);
    
    % adjust windows and keep indices
    for I = 1:nfilts
        tmpidx   = find(wts{I}>=lin);
        idx(I,:) = [min(tmpidx),max(tmpidx)];
        wts{I}   = wts{I}(tmpidx);
    end
end
