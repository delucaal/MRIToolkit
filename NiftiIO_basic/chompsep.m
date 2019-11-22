%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function str = chompsep(str)
%CHOMPSEP  Remove file separator at end of string.
%		STR = CHOMPSEP(STR) returns the string STR with the file separator at
%		the end of the string removed (if existing). 
%
%		Example:
%		str1 = chompseq('/usr/local/');
%		str2 = chompseq('C:\Program Files\');
%
%		Markus Buehren
%		Last modified 05.04.2009
%
%		See also CONCATPATH.

if ~isempty(str) && str(end) == filesep
  str(end) = '';
end