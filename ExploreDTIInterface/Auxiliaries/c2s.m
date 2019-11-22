%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%



% From ExploreDTI: converts cartesian coordinates to spherical
function S = c2s(C)
% Originally written from Ben Jeurissen under the supervision of Alexander
% Leemans
norm = sqrt(sum(C.^2, 2));
S(:,1) = acos(C(:,3)./norm);
S(:,2) = atan2(C(:,2), C(:,1));
end