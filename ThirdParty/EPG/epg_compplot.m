function epg_compplot(FZ,N)
%
%function epg_compplot(FZ,N)
%
%	Function plots the Mx, My and Mz components
%	across a voxel given the F and Z states.
%
%	FZ = FZ state matrix
%	N = number of states to plot (default = 50)
%

if (nargin < 2) N=50; end;
M = epg_FZ2spins(FZ,N);

z = 1:N;

subplot(1,3,1);
plot(M(1,:),z/N,'b-');
hold on; plot([0 0],[0 1],'k-'); hold off;
xlabel('Mx');
ylabel('Position in Voxel');
grid on;
title('Mx');


subplot(1,3,2);
plot(M(2,:),z/N,'r-');
hold on; plot([0 0],[0 1],'k-'); hold off;
xlabel('My');
ylabel('Position in Voxel');
grid on;
title('My');


subplot(1,3,3);
plot(M(3,:),z/N,'g-');
hold on; plot([0 0],[0 1],'k-'); hold off;
xlabel('Mz');
ylabel('Position in Voxel');
grid on;
title('Mz');



