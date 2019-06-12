%
%
%	Simulation of SPGR using EPG.  Yes, you can do this!!
%
%

Ntr = 600;	% Number of TRs
Nf = 144;
Nstates = 40;	% Number of states to simulate (easier to keep constant
		% if showing state evolution.

flipang = 20;	% Degrees
TR = 0.01;	% s
T1 = .5;	% s
T2 = 0.1;	% s
phaseinc=117;	% Degreees

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization.

Pstore = zeros(2*Nstates,Ntr);	% Matrix to store for image of evolution.
Zstore = zeros(Nstates,Ntr);	% Matrix to store for image of evolution.

% -- excitation

flipphase = 0;			% Start with flip angle along x (arbitrary)
flipphaseincr = 0;		% Phase increment 0 to start.
flipphaseincrincr = phaseinc;	% Phase increment-increment of 117 deg.

for k = 1:Ntr

  flipphase = flipphase + pi/180*(flipphaseincr);	% Incr. phase
  flipphaseincr = flipphaseincr + flipphaseincrincr;	% Increment increment!

  fp(k)=flipphase;
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

figure(2);
magphase(sd);

figure(3);
Np = 200;

for q=k-Nf:k
  Fp = [flipud(Pstore(1:Nstates,q)).'];
  Fn = [Pstore(Nstates:2*Nstates-1,q).']; Fn(1)=conj(Fn(1));
  Zn = [Zstore(:,q).'];
  M = epg_FZ2spins([Fp;Fn;Zn],Np)*Np;
  S = exp(-i*fp(q))*(M(1,:)+i*M(2,:));
  wmagphase(S);
  subplot(2,1,1);
  hold on;
  h = plot([1,200],[1 1]*exrecsignal(T1,T2,0,TR,flipang),'y--');
  plot([1,200],[1 1]*abs(mean(S)),'w--');
  hold off;
  setprops
  subplot(2,1,1);

  tt = sprintf('Signal Across Voxel, %g^\\circ Increment',phaseinc);
  h = title(tt);
  set(h,'FontSize',24);
  set(h,'Color',[1 1 1]*.8);
  drawnow;


  f = getframe(gcf);
  [fim,map] = frame2im(f);
  fn = sprintf('%s_%04d.tif','/Users/brian/tmp/ims/spgr',q+Nf-Ntr);
  imwrite(fim,fn,'TIFF');

end;

