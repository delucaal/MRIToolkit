%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



%	epg_bssfp
%
%	This is a simple example of how EPG can be
%	used to simulate a balanced SSFP profile.
%
%	By itself, EPG should not use this, as there is no
%	gradient spoiling in bSSFP.  However, the conversion
%	from EPG states to spin states allows the profile to
%	be plotted, and thus may allow "sub-states" to be used.

%	First, Create a gradient-spoiled sequence
%
P = [0 0 1]';		% Equilibrium
for k=1:100
  P = epg_grelax(P,1,.2,.01);	% Gradient spoiler.
  P = epg_grelax(P,1,.2,.01);	% Again...!  Show 2 cycles in the end .	%
  P = epg_rf(P,pi/6,pi/2);	% 30 degree tip.	% Sample signal here
end;

%	Convert to spins across voxel, which give the signal
%	vs off-resonance, essentially.

M = epg_FZ2spins(P);	% There are many redundant (0) states here that could	
			% be trimmed, but we actually want a high N to
			% show the effect of the simulation.

Msig = M(1,:)+M(2,:)*i;	% Signal
Msig = Msig*4*200;	% Show full signal.

subplot(2,1,1);
plot(abs(Msig));
title('Magnitude');
grid;
subplot(2,1,2);
plot(angle(Msig)/pi);
title('Phase / \pi');
grid;