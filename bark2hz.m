function hz=bark2hz(z)
%       BARK2HZ         Converts frequencies Bark to Hertz (Hz)
%	function hz=bark2hz(z)
%         Traunmueller-formula    for    z >  2 Bark
%         linear mapping          for    z <= 2 Bark
%
%	Author:   Kyrill, Oct. 1996
%                 Kyrill, March 1997   (linear mapping added)
%

% Hynek's formula (2003-04-11 dpwe@ee.columbia.edu)
% (taken from rasta/audspec.c)
hz = 600 * sinh(z/6);
