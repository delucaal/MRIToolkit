%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function toggle_elastix_options_in_parametersfile(filein,fileout,options)
% OPTIONS = CELL {'Key','NewValue',(Optional)'Write if not found' (1)}
tp = filein;
orig_file = fopen(tp,'rt');
new_file = fopen('tmp.txt','wt');
found_options = zeros(length(options),1);

while true
    line = fgetl(orig_file);
    if(feof(orig_file))
        break
    end
    nof_match = 0;
    for ij=1:length(found_options)
        if(~isempty(strfind(line,options{ij}{1})))
            if(found_options(ij) == 0)
                fprintf(new_file,'%s\n',options{ij}{2});
            end
            found_options(ij) = found_options(ij)+1;
            nof_match = nof_match+1;
            break
        end       
    end
    if(nof_match == 0)
        % NOTHING TO CHANGE, JUST COPY THE LINE
        fprintf(new_file,'%s\n',line);
    end

end

for ij=1:length(found_options)
    % CHECK IF SOME KEYS ARE MISSING AND SHOULD BE WRITTEN
    if(found_options(ij) > 0)
        continue
    end
    if(length(options{ij}) > 2 && options{ij}{3} == 1)
        fprintf(new_file,'%s\n',options{ij}{2});
        found_options(ij) = 1;
    end
end

fclose(orig_file);
fclose(new_file);

if(strcmp(fileout,filein))
    % IN CASE OF OVERWRITING SAVE OLD VERSION
    copyfile(filein, [filein '.old']);
end
copyfile('tmp.txt', fileout);
delete('tmp.txt');

end