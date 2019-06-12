function FZ = epg_gradanim(FZ,ns,nframes)
%
%	Animates the effect of gradient on F states
%	by plotting n intermediate frames between current
%	state and state that is achieved by one gradient pulse.
%
%	INPUT
%		FZ = F and Z state matrix.
%		ns = (optional) number of spins to show
%		nframes = (optional) number of intermediate states to show.
%
%	OUTPUT
%		FZ is 'advanced' to next state.


if (nargin < 2)
  ns = 50;	% Number of spins
end;

if ((nargin < 3) || (nframes<2))
  nframes= 20;
end;

for k=1:nframes
  epg_plot3DF(FZ,[],ns,(k-1)/(nframes-1));
  drawnow;
end;




