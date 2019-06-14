
%%% Distributed under the terms of CC BY-NC-ND (https://creativecommons.org/licenses) %%%
%
%function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)
%	Propagate EPG states through an RF rotation of 
%	alpha, with phase phi (both radians).
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		alpha = flip angle in radians.
%		phi = angle of rotation axis from Mx (radians).
%
%       OUTPUT:
%               FpFmZ = Updated FpFmZ state.
%		RR = RF rotation matrix (3x3).
%
%	SEE ALSO:
%		epg_grad, epg_grelax
%
%	B.Hargreaves.
%
function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)
% fprintf('epg_rf: %.2f - %.2f \n',alpha,phi);
% -- From Weigel at al, JMR 205(2010)276-285, Eq. 8.

if (abs(alpha)>2*pi) warning('epg_rf:  Flip angle should be in radians!'); end;

if (nargin < 3) warning('Rotation axis not specified - assuming -My'); 
	phi=-pi/2; 
end;

% -- Rotation matrix

RR = [(cos(alpha/2))^2 exp(2*i*phi)*(sin(alpha/2))^2 -i*exp(i*phi)*sin(alpha);
      exp(-2*i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 i*exp(-i*phi)*sin(alpha);
      -i/2*exp(-i*phi)*sin(alpha) i/2*exp(i*phi)*sin(alpha)      cos(alpha)];


FpFmZ = RR * FpFmZ;
