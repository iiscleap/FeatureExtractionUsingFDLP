function p = fdlpar(x,param)
%*****************************************************************
% USAGE P = FDLPAR(X,PARAM)
% FDLPAR Fit frequency-domain linear prediction model and return the poles
% Follows a Bark/Mel/Linear-Scale decomposition
% Function performs DCT of the signal, windows it into sub-bands and makes
% an autoregressive model of the sub-band envelopes
%*****************************************************************
% Sriram Ganapathy
% Center of Language and Speech Processing 
% Johns Hopkins University
% ganapathy@jhu.edu
%*****************************************************************
% 17-Jan-2012
% See the file COPYING for the licence associated with this software.
%*****************************************************************

if nargin < 2;  error ('NOT ENOUGH INPUT ARGUMENTS'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------- DCT and Sub-band Windowing ----------

if param.dct_high == 4000                   % Keep the DCT High to max value if this value is not set
    param.dct_high = param.fs/2;
end
if param.wiener 
   NIS = check_VAD(x,param.fs);
   NIS = find(NIS == 0);
   len_orig=length(x);	
end

lo_freq = floor(2*length(x)*param.dct_low/param.fs)+1 ; 
hi_freq = round(2*length(x)*param.dct_high/param.fs) ;

% Take the DCT
x = dct(x); 
% Restrict the DCT to the required region of interest
x = x(lo_freq:hi_freq);
% Get the frame length for FDLP input 
flen = size(x,1);

% Make the sub-band weights
if param.axis == 1              % Bark Scale of Decomposition     
    [wts,idx] = barkweights(flen,param.fs);
    
elseif  param.axis == 2         % Mel Scale of Decomposition
    [wts,idx] = melweights(flen,param.fs);
    
elseif param.axis == 3          % Linear Scale of Decomposition 
    nbands = 96;                % 96 sub-bands performs the best in reverb    
    [wts,idx]  = unif_rect_wind_fixed(nbands,flen);	      % Determine the sub-band DCT windows
else 
    error ('Unknown Frequency Decomposition');
end

% Define the model order 
fp = round(param.order*flen/param.fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------  Pole Estimation ----------------------

nb = size(idx,1);           % Number of sub-bands ...
p  = cell(nb,1);

% Time envelope estimation per band 
for I = 1+param.skip_bands:nb,
    % Apply the weights and get poles 
      tmpx = full(diag(sparse(wts{I}))*x(idx(I,1):idx(I,2),:));
    %******** Finally, FDLP ********% 
      if ~(param.wiener)
	      p{I} = do_lpc(tmpx,fp,param.gain_norm,param.lptype);  % FDLP !!!   
      else
	      p{I} = do_lpc_wiener2(tmpx,param.fs,len_orig, fp);  % FDLP !!!
      end	
end
