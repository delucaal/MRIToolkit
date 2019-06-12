%
%
%	Simulation of SPGR using EPG.  Yes, you can do this!!
%	(Version w/o phase)
%
%

Ntr = 200;	% Number of TRs
Nstates = 20;	% Number of states to simulate (easier to keep constant
		% if showing state evolution.

flipang = 20;	% Degrees
TR = 0.01;	% s
T1 = 1;		% s
T2 = 0.2;	% s

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization.

Pstore = zeros(2*Nstates,Ntr);	% Matrix to store for image of evolution.
Zstore = zeros(Nstates,Ntr);	% Matrix to store for image of evolution.

% -- excitation

flipphase = 0;			% Start with flip angle along x (arbitrary)
flipphaseincr = 0;		% Phase increment 0 to start.
flipphaseincrincr = 117;	% Phase increment-increment of 117 deg.

for k = 1:Ntr

  flipphase = flipphase + pi/180*(flipphaseincr);	% Incr. phase
  flipphaseincr = flipphaseincr + flipphaseincrincr;	% Increment increment!

  P = epg_rf(P,pi/180*flipang,flipphase);		% RF pulse
  s(k) = P(1,1);      			% Signal is F0 state.
  sd(k) = s(k)*exp(-i*flipphase);	% Phase-Demodulated signal.

  % -- Keep states for evolution diagram
  Pstore(Nstates:2*Nstates-1,k) = P(2,:).'; % Put in negative states
  Pstore(1:Nstates,k) = flipud(P(1,:).');  % Put in positive, overwrite center.
  Zstore(:,k) = P(3,:).';

  % -- Simulate relaxation and spoiler gradient
  P = epg_grelax(P,T1,T2,TR,1,0,1,1);   % spoiler gradient, relaxation.

end;

figure(1);
plotstate = cat(1,Pstore,Zstore);
dispim(plotstate,0,0.2);
xlabel('Time'); ylabel('F(top) and Z(bottom) States');


% Calculate theoretical perfectly-spoiled signal
s = sin(flipang*pi/180);
c = cos(flipang*pi/180);
TE=0;
E = exp(-TE/T2);
R = exp(-(TR)/T1);

if (abs(1-R*c)<.0001)
  sig = 1;			
else
  sig = (1-R)*E*s/(1-R*c);
end;


figure(2);
plot([1:length(sd)],abs(sd),'b-',[1:length(sd)],ones(1,length(sd))*abs(sig), 'r--');
xlabel('Repetition #');
ylabel('Signal');
legend('EPG Signal','Perfect Spoiling');
title('RF-Spoiled Evolution');
grid on;

