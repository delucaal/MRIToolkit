%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function dateStr = translatedatestr(dateStr)
%TRANSLATEDATESTR  Translate german date string to english version.
%		STR = TRANSLATEDATESTR(STR) converts a german date string like 
%		  13-M�-2006 15:55:00
%		to the english version
%		  13-Mar-2006 15:55:00.
%		This is needed on some systems if function DIR returns german date
%		strings.
%
%		Markus Buehren
%
%		See also DATENUM2.

dateStr = strrep(dateStr, 'Mrz', 'Mar');
dateStr = strrep(dateStr, 'M�', 'Mar');
dateStr = strrep(dateStr, 'Mai', 'May');
dateStr = strrep(dateStr, 'Okt', 'Oct');
dateStr = strrep(dateStr, 'Dez', 'Dec');