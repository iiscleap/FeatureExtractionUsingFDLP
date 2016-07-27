function a=do_lpc_wiener(signal,fs,flen,fp,NIS)
% ***************************************************************
% USAGE
% output=do_lpc_wiener(signal,fs,flen,fp,NIS)
% Wiener filtering using noise estimates obtained from the ETSI VAD
% Implementation adapted from Rainer Martin IEEE SP 2007
% Feb 2011
% ***************************************************************
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

W=fix(.025*fs);                % Window length is 25 ms
SP=0.4;                        % Shift percentage is 25ms - No overlap

wnd=rectwin(W);              
alpha=0.9; 

L = length(signal);
ENV_cmplx=fft(signal,2*flen-1);             % Full Temporal Envelope        
ENV=abs(ENV_cmplx(1:fix(end/2)+1,:)).^2;    %Full length temporal envelope

Y=segment(ENV,W,SP,wnd);                    % This function chops the signal into frames

numberOfFrames=size(Y,2);
Pn=mean(Y(:,NIS)')';                        % UDIM %initial Noise ENV variance 
X=zeros(size(Y));                           % Initialize X (memory allocation)
G=zeros(size(Y));                           % Initialize G (memory allocation)
zeta = zeros(size(Y));
gamma = zeros(size(Y));

for i=1:numberOfFrames
    if i > 1
        gamma(:,i) = Y(:,i)./Pn;
        zeta(:,i) = alpha*(X(:,i-1)./Pn) + (1-alpha) *  ( (gamma(:,i)-1) ) ;
    else
        gamma(:,i) = Y(:,i)./Pn; 	
        zeta(:,i) = (1-alpha) * ( (gamma(:,i)-1) );
    end
    G(:,i) = zeta(:,i)./(1+zeta(:,i));                 % Wiener Filter gain
 
    X(:,i)= G(:,i).^2.*Y(:,i);                        %Obtain the new Cleaned value
end

YPhase = zeros(size(X));
ENV_output=OverlapAdd2(X,YPhase,W,SP*W); %Overlap-add Synthesis of speech
ENV_output(length(ENV_output)+1:flen) = ENV(length(ENV_output)+1:flen); % Put the original envelope for the chopped frame

ENV_cmplx_op = ENV_output;
if mod(2*flen-1,2) %if FreqResol is odd
   ENV_cmplx_op =[ENV_cmplx_op;flipud(conj(ENV_cmplx_op(2:end,:)))];
else
    ENV_cmplx_op=[ENV_cmplx_op;flipud(conj(ENV_cmplx_op(2:end-1,:)))];
end
% '.' is just accelerated conjugation
R = real(ifft((ENV_cmplx_op/L)));

% Do it fast with levinson (works with complex)
[a,g2] = levinson(R,fp);

% Filter per column
a = a.';

% Get rid of the nasty imaginary roundoff if x is real
if isreal(signal)
    a = real(a);
end




function ReconstructedSignal=OverlapAdd2(XNEW,yphase,windowLen,ShiftLen);

%Y=OverlapAdd(X,A,W,S);
%Y is the signal reconstructed signal from its spectrogram. X is a matrix
%with each column being the fft of a segment of signal. A is the phase
%angle of the spectrum which should have the same dimension as X. if it is
%not given the phase angle of X is used which in the case of real values is
%zero (assuming that its the magnitude). W is the window length of time
%domain segments if not given the length is assumed to be twice as long as
%fft window length. S is the shift length of the segmentation process ( for
%example in the case of non overlapping signals it is equal to W and in the
%case of %50 overlap is equal to W/2. if not givven W/2 is used. Y is the
%reconstructed time domain signal.

if nargin<2
    yphase=angle(XNEW);
end
if nargin<3
    windowLen=size(XNEW,1)*2;
end
if nargin<4
    ShiftLen=windowLen/2;
end
if fix(ShiftLen)~=ShiftLen
    ShiftLen=fix(ShiftLen);
    disp('The shift length have to be an integer as it is the number of samples.')
    disp(['shift length is fixed to ' num2str(ShiftLen)])
end

[FreqRes FrameNum]=size(XNEW);

Spec=XNEW.*exp(j*yphase);
Win = repmat(rectwin(windowLen),1,FrameNum);
sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);
inv_win = zeros((FrameNum-1)*ShiftLen+windowLen,1);
weight=sig;
for i=1:FrameNum
    start=(i-1)*ShiftLen+1;    
    spec=Spec(:,i);
    
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+spec;   
    inv_win(start:start+windowLen-1)=inv_win(start:start+windowLen-1)+Win(:,i);
end
ReconstructedSignal=sig./inv_win;

function Seg=segment(signal,W,SP,Window)

% SEGMENT chops a signal to overlapping windowed segments
% A= SEGMENT(X,W,SP,WIN) returns a matrix which its columns are segmented
% and windowed frames of the input one dimentional signal, X. W is the
% number of samples per window, default value W=256. SP is the shift
% percentage, default value SP=0.4. WIN is the window that is multiplied by
% each segment and its length should be W. the default window is hamming
% window.

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

