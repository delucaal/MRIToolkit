%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of CC BY-NC-ND (https://creativecommons.org/licenses) %%%
function run_folder = get_executed_file_path() 
%     run_folder = mfilename('fullpath');
    run_folder = dbstack('-completenames');
    run_folder = run_folder(2).file;
    if(~isempty(run_folder))
        if(ispc)
            slash_location = strfind(run_folder,'\');
        else
            slash_location = strfind(run_folder,'/');
        end
        run_folder = run_folder(1:slash_location(end)-1);
    else
        run_folder = pwd;
    end 
end