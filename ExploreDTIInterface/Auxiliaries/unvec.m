
% From ExploreDTI: invert the previous operation and restore the original
% matrix shape
function y = unvec(x,mask)
%
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

dims = [size(mask,1) size(mask,2) size(mask,3)];
if isfloat(x)
    y = NaN([dims(1) dims(2) dims(3) size(x,1)], class(x));
else
    y = zeros([dims(1) dims(2) dims(3) size(x,1)], class(x));
end

for k = 1:size(x,1)
    if isfloat(x)
        Dummy = NaN(dims, class(x));
    else
        Dummy = zeros(dims, class(x));
    end
    Dummy(mask) = x(k,:);
    y(:,:,:,k) = Dummy;
end
end