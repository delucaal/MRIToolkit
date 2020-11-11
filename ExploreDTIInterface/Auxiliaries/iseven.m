function ret = iseven(x)
%
% Returns true if is is even (false otherwise).
%
% Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
%
ret = (mod(x,2) == 0);
end