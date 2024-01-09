%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

classdef ImageRegistrations < handle
    methods(Static)
        function elastix_cmd = BasicElastixCall()
            global MRIToolkit;
            if(isfield(MRIToolkit,'Elastix'))
                elastix_path = MRIToolkit.Elastix.Location;
                elastix_cmd = MRIToolkit.Elastix.ElastixCMD;
                if(~ispc)
                    if(~ismac)
                        elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' elastix_cmd];
                    else
                        elastix_cmd = ['DYLD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'elastix_Mac_64'];
                    end
                end
            end
        end

        function transformix_cmd = BasicTransformixCall()
            global MRIToolkit;
            if(isfield(MRIToolkit,'Elastix'))
                elastix_path = MRIToolkit.Elastix.Location;
                transformix_cmd = MRIToolkit.Elastix.TransformixCMD;
                if(~ispc)
                    if(~ismac)
                        transformix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' transformix_cmd];
                    else
                        transformix_cmd = ['DYLD_LIBRARY_PATH=' elastix_path ' ' transformix_cmd];
                    end
                end
            end
        end

        function TransfPar = ElastixRotationMatrix(elastix_transform_file)
            TransfPar = [];
            
            if(exist(elastix_transform_file,"file") < 1)
                disp('Cannot find the specified file')
                return
            end

            f = fopen(elastix_transform_file,'rt');
            
            parameters = {};
            while(feof(f) == false)
                l = fgetl(f);
                if(contains(l,'(TransformParameters'))
                    l = strrep(l,'(','');
                    l = strrep(l,')','');
                    parameters = strsplit(l);
                    parameters = parameters(2:end);
                    break
                end
            end
            
            fclose(f);
            
            TransfPar = [
                str2double(parameters{2})*(180/pi) str2double(parameters{1})*(180/pi) str2double(parameters{3})*(180/pi);
                str2double(parameters{11}) str2double(parameters{10}) str2double(parameters{12});
                str2double(parameters{8}) str2double(parameters{7}) str2double(parameters{9});
                str2double(parameters{5})  str2double(parameters{4}) str2double(parameters{6})
            ];
            TransfPar(isnan(TransfPar)) = 0;
        end

        function DownloadElastix()
            global MRIToolkit;
            download_url = '';
            if(ispc)
                download_url = 'https://github.com/SuperElastix/elastix/releases/download/5.1.0/elastix-5.1.0-win64.zip';
            elseif(ismac)
                download_url = 'https://github.com/SuperElastix/elastix/releases/download/5.1.0/elastix-5.1.0-mac.zip';
            elseif(isunix)
                download_url = 'https://github.com/SuperElastix/elastix/releases/download/5.1.0/elastix-5.1.0-linux.zip';
            end
            dest_folder = userpath;
            if(isempty(dest_folder))
                dest_folder = MRIToolkit.RootFolder;
            end
            urlwrite (download_url, fullfile(dest_folder,'Elastix.zip'));
            unzip(fullfile(dest_folder,'Elastix.zip'),fullfile(dest_folder,'Elastix'));
            delete(fullfile(dest_folder,'Elastix.zip'));
            if(ismac || isunix)
                system(['mv ' fullfile(dest_folder,'Elastix','bin','*') ' ' fullfile(dest_folder,'Elastix')]);
                system(['mv ' fullfile(dest_folder,'Elastix','lib','*') ' ' fullfile(dest_folder,'Elastix')]);
                system(['chmod +x ' fullfile(dest_folder,'Elastix','elastix')]);
                system(['chmod +x ' fullfile(dest_folder,'Elastix','transformix')]);
            end

            MRIToolkitGetSetVars('MRIToolkit.Elastix.Location',fullfile(dest_folder,'Elastix'));
            MRIToolkit.Elastix.Location = fullfile(dest_folder,'Elastix');
            MRIToolkit.Elastix.ElastixCMD = fullfile(MRIToolkit.Elastix.Location,'elastix');
            MRIToolkit.Elastix.TransformixCMD = fullfile(MRIToolkit.Elastix.Location,'transformix');            
        end

    end
end