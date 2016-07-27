function ceps = fdlp_feat(samples,config_file)
%*****************************************************************
% Function to extract FDLP based features
% Usage OUTPUT = FDLP_FEAT(SAMPLES,CONFIG-FILE)
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

if nargin < 1
    error ('Unknown Input - OUTPUT = FDLP_FEAT(SAMPLES,CONFIG-FILE');    
end

Time = clock ;
disp ( ' ****************************************** ' );
disp( ['FDLP Feature Extraction Invoked On ', num2str(datestr(now))]);

if nargin < 2
    disp ('Using default configuration from matlab.config file');
    config_file = 'matlab.config';
end
tic ;
%*****************************************************************
param = read_config_file(config_file) ;
%*****************************************************************

if max(abs(samples)) < 1
    samples = samples * 2^15;           % Making the samples to raw format
end

A = samples(:);
sr = param.fs;
param.flen= (param.fr_len/1000)*sr;           % frame length corresponding to 25ms
param.fhop= (param.fr_shift/1000)*sr;         % frame overlap corresponding to 10ms
fnum = floor((length(A)-param.flen)/param.fhop)+1;
send = (fnum-1)*param.fhop + param.flen;
A = A(1:send);
if param.fullsig
   fdlpwin = send;	
else
   fdlpwin = param.fdlplen*sr;             % FDLP window length
end

fdlpolap = ceil((param.fr_len-param.fr_shift)/10)*sr/100;
[X,add_samp] = frame_new(A,fdlpwin,fdlpolap);
ceps = [];

for i = 1 :size(X,2)                    % Go over each FDLP window
    x = X(:,i);
    % Now lets dither (make sure the original waves are not normalized!)
    x = ditherit(x);
    x = x - 0.97* [0 ; x(1:end-1)];                 % Pre-emphasis
    % FDLP processing starts here
    temp = do_feats_for_seg(x,param);
    ceps = [ceps temp];
end
ceps = ceps(:,1:fnum);
toc;
disp ( ['Completed ' num2str(fnum) ' feature vectors of dim ' num2str(size(ceps,1)) ...
    ' for a ' num2str(send/param.fs) ' sec file'])
disp ( ' ****************************************** ' );






