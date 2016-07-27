function param = read_config_file(config_file)
%*****************************************************************
% Function to parse the input config file and read the parameters.
% Config file must contain lines in format NAME = VALUE
% Detailed information about the parameters present in README
%*****************************************************************
% Sriram Ganapathy
% Center of Language and Speech Processing 
% Johns Hopkins University
% ganapathy@jhu.edu
%*****************************************************************
% 11-Jan-2012
% See the file COPYING for the licence associated with this software.


%*****************************************************************
%                   Default Parameters
%*****************************************************************
if nargin < 1 
    default = 1 ;
    warning ('No Config File Defined : Using Default Parameters');
else 
    default = 0;
end

% Signal Parameters
param.fs = 8000 ;             % Sampling Rate
param.fdlplen = 1 ;           % FDLP frame length in seconds (1 second)
param.fullsig = 0;	      % Using the entire signal in one frame length (useful for short speech files of 2-3s)	
% Type of Feature
param.type = 1;               % Type of features (fdlps=1/fdlpm=2/fdlp_plp2=3)
param.axis = 1;               % Frequency Axis   (bark=1/mel=2/linear=3)
param.mel_warp =0 ;           % Applies only for converting linear to mel-scale  
param.dct_low = 0;            % Begin of Frequency Axis in DCT domain
param.dct_high = 4000;        % End of Frequency Axis in DCT domain
param.skip_bands = 0 ;        % Skip the initial bands in feature computation

% Feature configuration
param.gain_norm = 0;          % Gain normalization 
param.order     = 160;        % Model order per sub-band per second
param.lptype    = 1 ;         % Type of linear prediction (autocorr=1/ls=2)
param.wiener    = 0 ;         % Apply Wiener filtering in feature

% Modulation feature configuration
param.num_temp_ceps = 14;     % Number of modulation components per band  
param.include_adapt = 1;      % Including static and adaptive compression
                              % Requires mex file adapt_m
                                                            
% Spectral feature parameters
param.num_spec_ceps   = 13;        % Number of cepstral components
param.include_c0      = 1;         % Flag set to 1 for using C1-C13
param.fr_len          = 25;        % Frame length for spectral frame (ms)
param.fr_shift        = 10;        % Frame shift for spectral fram (ms)   
param.delta           = 1;         % Apply Delta and Acceleration 
                                   % Only for Spectral Features    
%*****************************************************************

if default 
    return 
end

%*****************************************************************
%                   Config File Parameters
%*****************************************************************

[name,delim,value]=textread( config_file, '%s %s %d' );

N = length(name);

for i = 1 : N
    switch lower(name{i})
        case 'sample_rate'
            param.fs = value(i);
        case 'fdlplen'
            param.fdlplen = value(i);
	case 'flag_fullsiglen'
            param.fullsig = value(i);	
	case 'flag_feat_type'
            param.type = value(i);
        case 'flag_axis'
            param.axis = value(i);
        case 'flag_gain_norm'
            param.gain_norm = value(i);
        case 'model_order'
            param.model_order = value(i);
        case 'lp_type'
            param.lptype = value(i);
        case 'flag_wiener'
            param.wiener = value(i);
        case 'temp_ceps'
            param.num_temp_ceps = value(i);
        case 'spec_ceps'
            param.num_spec_ceps = value(i);
        case 'flag_adapt'
            param.include_adapt = value(i);
        case 'flag_c0'
            param.include_c0 = value(i);
        case 'spec_fr_len'
            param.fr_len = value(i);
        case 'spec_fr_shift'
            param.fr_shift = value(i);
        case 'model_order'
	    param.order = value(i);	
	case 'flag_delta'
            param.delta = value(i);
        case 'dct_low'
            param.dct_low = value(i);
        case 'dct_high'
            param.dct_high = value(i);
        case 'mel_warp'
            param.mel_warp = value(i);
        case 'skip_bands'
            param.skip_bands = value(i);
    end
end


            



