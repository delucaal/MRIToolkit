%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% From ExploreDTI: converts spherical coordinates to cartesian
function C = s2c(S)
% Originally written from Ben Jeurissen under the supervision of Alexander
% Leemans
C = [ cos(S(:,2)).*sin(S(:,1)), ...
    sin(S(:,2)).*sin(S(:,1)), ...
    cos(S(:,1)) ];
end