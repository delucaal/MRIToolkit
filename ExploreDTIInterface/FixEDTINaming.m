file2scan = 'EDTI.m';
adapted_file = 'EDTI_autofix.m';

f = fopen(file2scan,'rt');
fo = fopen(adapted_file,'wt');

while(~feof(f))
   line = fgetl(f);
   if(contains(line,'(') && contains(line,')') && ~contains(line,'function'))
      % likely a function
      parts = strsplit(line,' ');
      for pid=1:length(parts)
         if(contains(parts{pid},'('))
            idx = strfind(parts{pid},'(');
            subpart = parts{pid}(1:idx-1);
            check_function = which(subpart);
            if(~isempty(check_function))
                if(contains(check_function,'ExploreDTI') && ~contains(check_function,'EDTI') && ...
                        contains(check_function,'.m') && ~contains(subpart,'.') && ~contains(subpart,'mex') ...
                        && ~contains(check_function,'MRIToolkit'))
%                     disp(['Possible EDTI function ' subpart ' in ' line]);
                    line = strrep(line,subpart,['EDTI_Library.' subpart]);
                    disp(['New line ' line]);
                end
            end
         end
      end
   end
   fprintf(fo,'%s%s',line,newline);
   
end

fclose(f);
fclose(fo);
