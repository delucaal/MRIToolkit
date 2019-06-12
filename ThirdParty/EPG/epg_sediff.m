%
%	Very simple simulation of Spin-echo diffusion.
%

D = 2*10^(-9);	% Diffusion (m^2/s)
gamma = 4258*2*pi;	% Gamma, Rad/s/G.

N = 200;	% Sequence repetitions.
TR = .01;	% 10ms


G = 8;		% G/m (Smallest gradient amplitude)
Tg = .02;	% Gradient Duration (s)
TR = .02;	% s - space between every RF and each gradient.
		% Note TE = 2* (2*TR+Tg)
T1 = 0.5;	% s
T2 = .038;	% s

ng = 50;

dk = gamma*G*Tg;
clear('S');

% ====== START SEQUENCE =========

for n = 1:ng	% Loop with successively bigger gradients.

  dkn = dk*n;	% Gradient "dk" for this scaling.
 
  F = [0;0;1];	% Equilibrium Magnetization.
		% [F+; F-; Z],  all longitudinal in Z0 state.

  F = epg_rf(F,pi/2,0);			% 90 RF Rotation.
  F = epg_grelax(F,T1,T2,TR,dkn,D,0,0);	% Relaxation, Diffusion, "TR"
  F = epg_grelax(F,T1,T2,Tg,dkn,D,1,0);	% Relaxation, Diffusion, Gradient
  F = epg_grelax(F,T1,T2,TR,dkn,D,0,0);	% Relaxation, Diffusion, "TR"
  
%   F = epg_rf(F,pi,pi/2);			% 180 RF Rotation.
  F = epg_rf(F,pi,0);			% 180 RF Rotation.

  F = epg_grelax(F,T1,T2,TR,dkn,D,0,0);	% Relaxation, Diffusion, "TR"
  F = epg_grelax(F,T1,T2,Tg,dkn,D,1,0);	% Relaxation, Diffusion, Gradient
  F = epg_grelax(F,T1,T2,TR,dkn,D,0,0);	% Relaxation, Diffusion, "TR"

  S(n) = F(1,1);			% Signal is F+(1) state.  

  % -- Calculate B value
  b_val(n) = ((gamma*n*G*Tg)^2)*((2*TR+Tg)-Tg/3);

end;

% -- Plot stuff.
db = b_val-b_val(1);

S = S./abs(S(1));
figure(1);
hold off;
plot(b_val/1000000,abs(S),'b-');	% Plot simulated signal
hold on;
a = axis; axis([a(1:2) 0 1]);
plot(db/1000000,exp(-db*D),'r--');	% Plot simple exp(-bD) prediction.
xlabel('B-Value (s/mm^2)');
ylabel('Signal');
legend('EPG Calculated','exp(-bD) Calculation');
title('Diffusion-Weighted Imaging');
setprops;
grid



 
