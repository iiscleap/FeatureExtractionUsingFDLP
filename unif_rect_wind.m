function [wts,idx]  = unif_rect_wind(wlen,whop,flen)
% -----------------------------------------------------
% Usage: [wts,idx]  = unif_rect_wind(wlen,whop,flen)
% function to obtain uniform rect windows ...
% of a given window size wlen and window hop whop ...
% Returns the sub-band window with square shape in wts
% The index of the sub-band DCTs is present in idx
% --------------------------------------------------

w = ones(1,flen);			% Sub-band window with square shape
[wind,addSamp] = frame_new(w,wlen,wlen-whop);	% Ibtain the sub-band windows
N_fr = size(wind,2) ; 
wts = cell(1,N_fr);
idx = zeros(N_fr,2);			% Array containing the indices for each sub-band with w.r.t full band DCT.
idx(1,:) = [1 wlen];			

for I  = 1 : N_fr			% Finding the starting and end-points of each sub-band.
    wts{I} = wind(:,I);
    if I > 1 
        idx(I,:) = idx(I-1,:)+whop;
    end
end


    
    
