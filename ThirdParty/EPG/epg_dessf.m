%
%function [S1,S2,F] = epg_dessf(T1,T2,TR,TE,alpha,G,Tg,D,N)
%
%	EPG Simulation of DESS (Double-Echo in Steady State)
%
%	Sequence is alpha - TE - Tx - Tg - Tx - TE  (Duration TR)
%
%	TE = echo time, Tg = gradient time, Tx = (TR-Tg)/2-TE (all seconds)
%	G = G/m, alpha=degrees.
%	T1,T2 = seconds
%	N = number of TRs, or fractional change in S1 and S2 at which to stop.
%	D = ADC, m^2/s
%
function [S1,S2,F] = epg_dessf(T1,T2,TR,TE,alpha,G,Tg,D,N)

if (nargin <9) nprec = 0.01; N=100; end;
if (N>1) 
  nprec = -1; % Don't stop!
else
  nprec=N; N=100;
end;	


alpha = pi/180*alpha;	% to radians.
gamma = 4258*2*pi;	% Gamma, Rad/G.


dk = gamma*G*Tg;


  F = [0;0;1];	% Equilibrium Magnetization.
		% [F+; F-; Z],  all longitudinal in Z0 state.

  F(1,10)=0;	% Allocate 10 states.
  noadd = 1;	% Don't add higher-order states.

  S1last =0;
  S2last =0;
  q = 1;
  done=0;
  while (done==0)
    F = epg_rf(F,alpha,0);			% 90 RF Rotation.

    F = epg_grelax(F,T1,T2,TE,dk ,D,0,noadd);	% Relaxation, Diffusion, to TE1
    SS1= F(1,1);			% Signal is F+(1) state.  

    F = epg_grelax(F,T1,T2,(TR-Tg)/2-TE,dk ,D,0,noadd);	% to start of G

    F = epg_grelax(F,T1,T2,Tg,dk ,D,1,noadd);	% Relaxation, Diffusion, Grad 1
  
    F = epg_grelax(F,T1,T2,(TR-Tg)/2-TE,dk ,D,0,noadd);	% end G to TE2
    SS2= F(1,1);			% Signal is F+(1) state.  
    F = epg_grelax(F,T1,T2,TE,dk ,D,0,noadd);	% Relax, Diffusion, TE2 to TR

    %F = epg_trim(F,0.0001);	% Remove higher-order states if < 0.0001
 
    % == Stopping criteria ==  
    q=q+1;
    if (q==N) 
      done=1; 
    else
      if (abs((abs(SS1)-abs(S1last))/abs(SS1)) < nprec )
        if (abs((abs(SS2)-abs(S2last))/abs(SS2)) < nprec )
	  done=1;
          tt=sprintf('DESS simulation converged in %d TRs',q); disp(tt);
        end;
      end;
    end;
    S1last=SS1;
    S2last=SS2;

  end;
  S1 = SS1;
  S2 = SS2;



 
