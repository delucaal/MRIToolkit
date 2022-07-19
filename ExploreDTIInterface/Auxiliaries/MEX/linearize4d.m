%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function y = linearize4d(x,mask)
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

if iscell(x)
    dims = length(x);
    maskSize = sum(mask(:));
    y = zeros([dims maskSize], class(x{1}));
    for k = 1:dims;
        y(k,:) = x{k}(mask(:));
    end   
else
    dims = size(x);
    maskSize = sum(mask(:));
    y = zeros([dims(4) maskSize], class(x));
    for k = 1:dims(4);
        Dummy = x(:,:,:,k);
        y(k,:) = Dummy(mask(:));
    end
end