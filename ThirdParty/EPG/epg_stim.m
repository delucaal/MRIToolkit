%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
%function [S,P] = epg_stim(flips)
%
%	Function uses EPG simulation to quickly calculate
%	the signal from a stimulated echo sequence defined
%	as
%		flip - gradient - flip - gradient ... flip - gradient
%
%	Relaxation is neglected, and flip angles (vector) are in degrees
%
%  	Example Try:  fplot('abs(epg_stim([x x x]))',[0,120]);
%


  
function [S,P] = epg_stim_calc(flips)
P = [0 0 1]';				% Equilibrium
for k=1:length(flips)		
  P = epg_rf(P,flips(k)*pi/180,pi/2);	% RF pulse
  P = epg_grelax(P,1,.2,0,1,0,1);	% Gradient
end;
S = P(1,1);				% Return signal from F0

