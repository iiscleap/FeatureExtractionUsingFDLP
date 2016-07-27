function [wts,idx]  = unif_rect_wind_fixed(nb,flen)

% Usage [wts,idx]  = unif_rect_wind_fixed(nb,flen)
% function to obtain uniform rect windows ...
% of fixed number of bands and overlap ...

% Set of fixed parameters 
whop = round(flen/(nb+3.5));		%  Determine the DCT window hop length
wlen = round(2.5*whop);			%  Determine the DCT window length	
[wts,idx]  = unif_rect_wind(wlen,whop,flen);

% Adjusting the last band length to adjust to the frame length.
if idx(end,end) > flen
    idx(end-1,end) = flen;
    wts{end-1} = ones(idx(end-1,2)-idx(end-1,1)+1,1);
    idx(end,:) = [];
    wts(end) = [];
elseif idx(end,end) < flen
    idx(end,end) = flen;
    wts{end} = ones(idx(end,2)-idx(end,1)+1,1);
end



    
    
