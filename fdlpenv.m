function env = fdlpenv(p,npts)
%*****************************************************************
% FDLPENV Calculate fdlp envelopes from FDLP poles
% Usage : env = fdlpenv(p,npts)
% input is one dimensional cell array in p
% npts is the number of points required in the envelope
% env is the output LP envelope
% Computed from LP polynomial coefficients and finding the inverse spectrum.  
%*****************************************************************
% Sriram Ganapathy
% Center of Language and Speech Processing 
% Johns Hopkins University
% ganapathy@jhu.edu
%*****************************************************************
% 25-Jan-2012
% See the file COPYING for the licence associated with this software.
%*****************************************************************

if nargin < 2
    npts = 1000;
end

% Find the size of the input signal
[na,nf] = size(p{1});

% How many fft points to calculate
nfft = 2*(max(npts,na)-1);

% Calculate the LP frequency response
h = fft(p{1},nfft);

% Use the positive half of the spectrum
h = h(1:(nfft/2)+1,:);

% Get the inverse spectrum
h = 1./h;

% Power spectrum 
h = (h.*conj(h));

% Correcting energy factor
env = 2*h;

