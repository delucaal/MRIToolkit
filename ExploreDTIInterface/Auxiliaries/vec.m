% From ExploreDTI: just cast a 4D matrix to a vector (without zeros)
function [y, mask] = vec(x,mask)
%
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

if ~exist('mask','var')
    mask = ~isnan(x(:,:,:,1));
end
y = zeros([size(x,4) sum(mask(:))], class(x));
for k = 1:size(x,4)
    Dummy = x(:,:,:,k);
    y(k,:) = Dummy(mask(:));
end
end