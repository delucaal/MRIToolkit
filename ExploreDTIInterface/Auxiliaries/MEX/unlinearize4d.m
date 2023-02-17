%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
function y = unlinearize4d(x,mask)
dims = size(mask);
y = NaN([dims(1) dims(2) dims(3) size(x,1)], class(x));
for k = 1:size(x,1)
    Dummy = NaN(size(mask), class(x));
    Dummy(mask) = x(k,:);
    y(:,:,:,k) = Dummy;
end
end