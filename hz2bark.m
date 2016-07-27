function z= hz2bark(f)
%       HZ2BARK         Converts frequencies Hertz (Hz) to Bark

% Inverse of Hynek's formula (see bark2hz)
% Written by Dan Ellis
z = 6 * asinh(f/600);

