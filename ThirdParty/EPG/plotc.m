%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% plotc(x,y) 
%
% plots the complex vector y as a function of x.
% plots real (blue --) ,imag (red :)  and magnitude (yellow -).
% returns the row as a vector. 

function plotc(x,y); 

if (nargin <2)
	y = x;
	x = 1:max(size(y));
	% color... plot(x,real(y),'b--',x,imag(y),'r:',x,abs(y),'y-');
	plot(x,real(y),'k--',x,imag(y),'k:',x,abs(y),'k-');
else
	% color... plot(x,real(y),'b--',x,imag(y),'r:',x,abs(y),'y-');
	plot(x,real(y),'k--',x,imag(y),'k:',x,abs(y),'k-');
end;