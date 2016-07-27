function flag_VAD = check_VAD(x,sr)

% FUnction to perform VAD similar to ETSI feature extraction
% Samples should be read from raw format file
% Details in ETSI ES 202 050 Document

% CONSTANTS

if sr ~= 8000
    x = resample(x,8000,sr); % Resample the test data to 8kHz for determining VAD information
end

sr = 8000;

NB_FRAME_THRESHOLD_LTE = 10;
LAMBDA_LTE = 0.97;
M =80;
SNR_THRESHOLD_UPD_LTE = 20;
ENERGY_FLOOR          = 80;
MIN_FRAME             = 10;

lambdaLTEhigherE      =0.99;
SNR_THRESHOLD_VAD     =15 ;

MIN_SPEECH_FRAME_HANGOVER =4;
HANGOVER             =15;

% Initialization
nbSpeechFrame =0;
meanEn =0;
hangOver =0;

% Frame the signal into 25ms windows with a shift of 10ms

flen = 0.025*sr;
W = hanning(flen);
SP =0.4;
x_fr = segment(x,flen,SP,W);
N_fr = size(x_fr,2);
flag_VAD = zeros(1,N_fr);

for t = 1 : N_fr
    x_cur = x_fr(:,t);
    if t < NB_FRAME_THRESHOLD_LTE
        lambdaLTE = 1 -1/t;
    else
         lambdaLTE = LAMBDA_LTE;
    end
    frameEn = 0.5 + 16/(log(2)) * (log((64 + sum(x_cur(end-79:end).^2))/64));
    
    if (frameEn - meanEn) < SNR_THRESHOLD_UPD_LTE || t < MIN_FRAME
        if frameEn < meanEn ||  t < MIN_FRAME
            meanEn = meanEn + (1-lambdaLTE)*(frameEn- meanEn);
        else
            meanEn = meanEn + (1-lambdaLTEhigherE)*(frameEn- meanEn);
        end
        if meanEn < ENERGY_FLOOR
            meanEn = ENERGY_FLOOR;
        end
    end
    
    if t > 4
        if (frameEn -meanEn) > SNR_THRESHOLD_VAD
            flag_VAD(t) = 1;
            nbSpeechFrame = nbSpeechFrame+1;
        else
            if nbSpeechFrame > MIN_SPEECH_FRAME_HANGOVER
                hangOver = HANGOVER;
            end
            nbSpeechFrame = 0;
            if hangOver ~= 0
                hangOver = hangOver - 1;
                flag_VAD(t) = 1;
            else
                flag_VAD(t) = 0;
            end
        end
    end
end          
    


function Seg=segment(signal,W,SP,Window)

% SEGMENT chops a signal to overlapping windowed segments
% A= SEGMENT(X,W,SP,WIN) returns a matrix which its columns are segmented
% and windowed frames of the input one dimentional signal, X. W is the
% number of samples per window, default value W=256. SP is the shift
% percentage, default value SP=0.4. WIN is the window that is multiplied by
% each segment and its length should be W. the default window is hamming
% window.
% 06-Sep-04
% Esfandiar Zavarehei

if nargin<3
    SP=.4;
end
if nargin<2
    W=256;
end
if nargin<4
    Window=hamming(W);
end
Window=Window(:); %make it a column vector

L=length(signal);
SP=fix(W.*SP);
N=fix((L-W)/SP +1); %number of segments

Index=(repmat(1:W,N,1)+repmat((0:(N-1))'*SP,1,W))';
hw=repmat(Window,1,N);
Seg=signal(Index).*hw;