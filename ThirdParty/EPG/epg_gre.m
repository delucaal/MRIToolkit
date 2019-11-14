%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
%
%	Very simple simulation of Gradient-Spoiled 
%	sequence using EPG methods.
%

N = 200;	% Number of Sequence repetitions.
TR = .01;	% 10ms

F = [0;0;1];	% Equilibrium Magnetization.
		% [F+; F-; Z],  all longitudinal in Z0 state.


S1 = zeros(N,1);
S2 = zeros(N,1);

for k=1:N
  F = epg_rf(F,pi/6,0);		% RF Rotation.
  S1(k) = F(1,1);		% Record "Signal"
  F = epg_grelax(F,1,.1,TR,1,0,1,0);	% Relaxation, spoiling
  S2(k) = F(1,1);		% Record "Signal"
end;
 
%	epg_plot - simple plot of F,Z states.
%	epg_plotcomp - plot of Mz,My,Mz components across a voxel
epg_plot(F);
plot([abs(S1) abs(S2)]);


 