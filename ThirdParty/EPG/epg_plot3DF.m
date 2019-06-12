%function epg_plot3DF(FZ,n,nm,frac)
%
%	Function plots the F state(s) as a 3D rendition
%	of spins within a voxel in the mx,my,mz space.
%
%	INPUT:
%		FZ = F/Z states, 3xM (Z are not used)
%		n = (optional) F state number to plot (0 to M-1), omit for all.
%		nm = (optional) number of spins in voxel to plot, omit for 30
%		frac = (optional) fraction of gradient precession before plot
%
function epg_plot3DF(FZ,n,nm,frac)

if ((nargin > 1) && (length(n)>0))
  if (n > size(FZ,2)) error('State must be within state matrix'); end;
  FZp = 0*FZ; FZp(:,n+1)=FZ(:,n+1); FZ = FZp;
end;

if (nargin < 3) nm = 30; end;
if (nargin < 4) frac = 0; end;



z = [1:nm]/nm;
M = epg_FZ2spins(FZ,nm,frac);

mx = max(abs(M(1,:)+i*M(2,:)));	% Maximum vector length
plot3(mx*[-1 1],[0 0],[0 0],'k-'); hold on;
plot3([0 0],mx*[-1 1],[0 0],'k-'); 
plot3([0 0],[0 0],[-0.1 1.1],'k-');

colormap('default');
c = colormap;
sz = size(c);
while(sz(1)<length(z))
  c = [c c];
  sz = size(c);
end;

for k=1:length(z)
  h=plot3([0 M(1,k)],[0 M(2,k)],z(k)*[1 1]);
  set(h,'Color',c(ceil(k/length(z)*sz(1)),:));
end;
hold off;

xlabel('Mx');
ylabel('My');
zlabel('Mz');









