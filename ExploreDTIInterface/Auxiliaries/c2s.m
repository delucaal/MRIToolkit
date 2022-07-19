%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% From ExploreDTI: converts cartesian coordinates to spherical
function S = c2s(C)
% Originally written from Ben Jeurissen under the supervision of Alexander
% Leemans
norm = sqrt(sum(C.^2, 2));
S(:,1) = acos(C(:,3)./norm);
S(:,2) = atan2(C(:,2), C(:,1));
end