function feats = do_feats_for_seg(x,param)
%*****************************************************************
% Function to extract FDLP based features for the current segment
% USAGE :  FEATS = DO_FEATS_FOR_SEG(X,PARAM)
% Cell array param defined using config file
% Config file must contain lines in format NAME = VALUE
% Detailed information about the parameters present in README
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
pad_fr = 2;
padlen=pad_fr*param.fhop;
% Padding the signal with 10 ms on both sides
x = [flipud(x(1:padlen)); x ; flipud(x(end -padlen+1:end))];

fnum = floor((length(x) - param.flen)/param.fhop) + 1 ;
send = (fnum-1)*param.fhop + param.flen ;
factor = 20*(param.fs/8000);           % Downsampling on the envelopes
fdlplen = floor(send/factor);

% FDLP pole extraction
p = fdlpar(x,param);
nb = size(p,1);                     % Number of Sub-bands

%*****************************************************************
%                   Envelope Generation
%*****************************************************************

fdlp_spec =zeros(nb,fdlplen);
for J = 1+param.skip_bands:(nb)
    % Envelope generation using polynomial interpolation
    ENV  =  fdlpenv(p(J),fdlplen);
    fdlp_spec(J,:) = ENV(:)';
end
clear ENV ;

if (param.type == 1 ) || (param.type == 3 ) 
    
    %*****************************************************************
    %                  I. Short-term Feature Extraction
    %*****************************************************************
    energy_bands=zeros((nb - param.skip_bands),fnum);
    wind = hamming(floor(param.flen/factor))';
    opt = 'nodelay';
    % Energy Integration
    for band = 1+param.skip_bands : size(fdlp_spec,1)
        [band_data,z,opt]=buffer(fdlp_spec(band,:),floor(param.flen/factor),floor((param.flen-param.fhop)/factor),opt);
        opt = 'nodelay';
        energy_bands((band-param.skip_bands),:) = (wind*band_data);
    end
    % Mel-warping the axis if needed
    if param.mel_warp
        energy_bands_new = audspec(energy_bands, sr,37,'mel');   % Convert fr. axis to mel scale.
        postspectrum = postaud(energy_bands_new,sr/2,'mel');
    else
        postspectrum = energy_bands;
    end
    
    % Cepstral Transform
    if param.type == 3
       % plp2 smoothing
	modelorder=min(18,size(energy_bands,1)); 	% TDLP Model Order
	% do LPC
	lpcas = dolpc(energy_bands, modelorder);	% TDLP
	% convert lpc to spectra 
	cepstra = lpc2ceps(lpcas, param.num_spec_ceps); % Convert to cepstra
        % One can also use lpc2spec to get the spectrogram back at this stage
 
    else
       if ~(param.include_c0)
        cepstra = spec2cep(postspectrum,param.num_spec_ceps+1);
        cepstra = cepstra(2:end,:);
       else
        cepstra = spec2cep(postspectrum,param.num_spec_ceps);
       end
    end
 
    cepstra = lifter(cepstra, 0.6);
    
    % Delta Double Delta ...
    if param.delta
        % Append deltas and double-deltas onto the cepstral vectors
        del = deltas(cepstra);
        
        % Double deltas are deltas applied twice with a shorter window
        ddel = deltas(deltas(cepstra,5),5);
        
        % Composite, 39-element feature vector, just like we use for speech recognition
        feats = [cepstra;del;ddel];
    else
        feats = cepstra;
    end
    
elseif param.type == 2
    %*****************************************************************
    %                II.  Modulation Feature Extraction
    %*****************************************************************
    trap = 10;                  % 10 FRAME context duration
    mirr_len = round(trap*param.fhop/factor);
    
    fdlpwin = floor((0.225*param.fs)/factor);         % Modulation spectrum Computation Window.
    fdlpolap = fdlpwin - round(param.fhop/factor);
    ceps_time = cell(1,nb);
    del = cell(1,nb);

    % Go over the bands and compute modulation spectra    
    for J = 1+param.skip_bands:(nb)
        ENV  =  fdlp_spec(J,:);
        ENV_log = log(ENV);
        ENV_log_new = [fliplr(ENV_log(1:mirr_len))  ENV_log fliplr(ENV_log(end-mirr_len+1:end))]; 
        ENV_log_fr = frame_new(ENV_log_new,fdlpwin,fdlpolap);
        ceps_time{J} = spec2cep_log(ENV_log_fr,param.num_temp_ceps);
        
        % Compute delta features using the Kollmeier Adaptative Compression model
        
        if param.include_adapt
            ENV_new = [repmat(ENV(:,1),1,floor(1000/factor)) ENV];    % Initial mirroring to initialize the adaptation model
            ENV_adpt = adapt_m(ENV_new/max(ENV_new),floor(param.fs/factor));
            ENV_adpt = ENV_adpt(round(end-send/factor+1):end);
            ENV_adpt_new = [fliplr(ENV_adpt(1:mirr_len))  ENV_adpt fliplr(ENV_adpt(end-mirr_len+1:end))]; 
            ENV_adpt_fr = frame_new(ENV_adpt_new,fdlpwin,fdlpolap);
            del{J} = spec2cep_log(ENV_adpt_fr,param.num_temp_ceps);                % Envelope already in the compressed domain
        end
    end
    feats = [];
    
    if param.include_adapt
        
        for J = 1+param.skip_bands : nb                             % get last 19 bands only
            feats = [feats ceps_time{J}' del{J}']; 
        end
    else
        for J = 1+param.skip_bands : nb
            feats = [feats ceps_time{J}']; 
        end
    end
    
    feats = feats'; % Inorder to have the notion of feature vectors
else
    error ('Unknown Feature Type, Choose 1-FDLPS or 2-FDLPM');
end

%%% Unpadding
feats = feats(:,pad_fr+1:end-pad_fr);
