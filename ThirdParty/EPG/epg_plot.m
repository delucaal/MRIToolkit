%
%	function epg_plot(FZ)
%
%	Function makes a simple plot of the F and Z states
%	in the 3xN array.
%
%	INPUT:
%		FZ	- 3xN matrix of Fp, Fm and Z states.
%
 
function epg_plot(FZ)


[y,x] = size(FZ);
xa = 0:x-1;

subplot(2,1,1);
plot(xa,abs(FZ(1,:)),'r-', xa,abs(FZ(2,:)),'b--', xa,abs(FZ(3,:)),'k-.');
legend('F+','F-','Z');
xlabel('State number');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(xa,angle(FZ(1,:)),'r-', xa,angle(FZ(2,:)),'b--', xa,angle(FZ(3,:)),'k-.');
legend('F+','F-','Z');
xlabel('State number');
ylabel('Phase');
grid on;
if (exist('setprops')) setprops; end;


