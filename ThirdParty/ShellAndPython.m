%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% This class implements methods to execute command line tools on different
% systems (BASH, ZSH, Powershell, WSL) and interface MRIToolkit with Python
% commands.

classdef ShellAndPython < handle
    methods(Static)
        function output = ExecuteSH(command,windows_wsl)
            if(ismac)
                cmd = ['. ~/.zshrc;. ~/.zprofile; ' command];
                [~,output] = system(cmd);
            elseif(isunix)
                cmd = ['. ~/.bashrc;. ~/.bash_profile; ' command];
                [~,output] = system(cmd);
            elseif(ispc)
                if(nargin < 2 || windows_wsl == false)
                    cmd = command;
                else
                    % To be appropriately fixed
                    cmd = ['. ~/.bashrc;. ~/.bash_profile; ' command];
                end
                [~,output] = system(cmd);
            end
        end

        function out = ExecuteCommandInCondaEnv(env,cmd)
            command = ['conda run -n ' env ' ' cmd];
            out = ShellAndPython.ExecuteSH(command);
        end

        function conda_loc = LocateAnaconda()
            global MRIToolkit;
            if(exist(MRIToolkit.Miniconda3Path,'dir') < 1)
                % Try to guess it
                if(ismac || isunix)
                    cmd = 'conda activate; which python';
                    output = ShellAndPython.ExecuteSH(cmd);
                    fs_hits = strfind(output,'/');
                    conda_loc = output(1:fs_hits(end));
                    MRIToolkit.Miniconda3Path = cond_loc;
                elseif(ispc)
                    disp('Cannot guess the path of Anaconda on your system. Please update your MRIToolkitDefiniteLocalVars.m file accordingly.')
                end
            else
                conda_loc = MRIToolkit.Miniconda3Path;
            end
        end

        function out = SetupScilpy()
            % Install scilpy if not already available
            cmd = 'conda run -n MRIToolkit which scil_compute_local_tracking.py';
            out = ShellAndPython.ExecuteSH(cmd);
            if(contains(out,'not found') || contains(out,'ERROR'))
                % Download and install
                disp('I will now attempt to install scilpy. If this fails, it might be due to missing prerequisites, please see also https://github.com/frheault/scilpy');
                cmd = 'conda activate MRIToolkit; pip install "packaging>=19.0"; pip install "numpy==1.23.*"; pip install "Cython==0.29.*"; pip install pyopencl';
                out = ShellAndPython.ExecuteSH(cmd);
                %if(contains(out,'Successfully installed Cython'))              
                cmd = ['cd ' tempdir2 ';rm -rf scilpy;git clone https://github.com/frheault/scilpy.git'];
                out = ShellAndPython.ExecuteSH(cmd);
                cmd = ['cd ' fullfile(tempdir2,'scilpy') ';conda run -n MRIToolkit pip install -e .'];
                out = ShellAndPython.ExecuteSH(cmd);
                if(contains(out,'ERROR'))
                    disp(out);
                else
                    disp('Success installing Scilpy');
                end
            end
        end

        function SetupANTs()
            global MRIToolkit;
            download_mac = 'https://github.com/ANTsX/ANTs/releases/download/v2.5.0/ants-2.5.0-macos-12-X64-clang.zip';
            download_linux = 'https://github.com/ANTsX/ANTs/releases/download/v2.5.0/ants-2.5.0-ubuntu-20.04-X64-gcc.zip';
            download_pc = 'https://github.com/ANTsX/ANTs/releases/download/v2.5.0/ants-2.5.0-windows-2022-X64-VS2019.zip';
            download_folder = fullfile(MRIToolkit.RootFolder,'ThirdParty','ANTs');
            if(exist(download_folder,'dir') < 1 || exist(fullfile(download_folder,'ants-2.5.0'),'dir') < 1)
                mkdir(download_folder);
                if(ispc)
                    urlwrite(download_pc,fullfile(download_folder,'ants.zip'));
                elseif(ismac)
                    urlwrite(download_mac,fullfile(download_folder,'ants.zip'));
                elseif(isunix)
                    urlwrite(download_linux,fullfile(download_folder,'ants.zip'));
                end
                unzip(fullfile(download_folder,'ants.zip'),download_folder);
                delete(fullfile(download_folder,'ants.zip'));
            else
                setenv('ANTSPATH',fullfile(download_folder,'ants-2.5.0','bin'));
                setenv('PATH',[getenv('PATH') ':' getenv('ANTSPATH')]);
            end
        end

        function SetupRecoBundlesX()
            global MRIToolkit;
            ShellAndPython.SetupANTs(); % Required
            download_link_atlas = 'https://zenodo.org/records/10103446/files/atlas.zip?download=1';
            download_link_config = 'https://zenodo.org/records/10103446/files/config.zip?download=1';
            download_folder = fullfile(MRIToolkit.RootFolder,'ThirdParty','RecoBundlesX');
            if(exist(download_folder,'dir') < 1 || isempty(dir(fullfile(download_folder,'*.json'))))
                try
                    mkdir(download_folder);
                catch
                end
                urlwrite(download_link_atlas,fullfile(download_folder,'atlas.zip'));
                urlwrite(download_link_config,fullfile(download_folder,'config.zip'));
                unzip(fullfile(download_folder,'atlas.zip'),download_folder);
                unzip(fullfile(download_folder,'config.zip'),download_folder);
                delete(fullfile(download_folder,'atlas.zip'));
                delete(fullfile(download_folder,'config.zip'));
            end
        end

        function out = DeleteDefaultPythonEnv()
            cmd = 'conda env remove -n MRIToolkit -y';
            out = ShellAndPython.ExecuteSH(cmd);            
        end

        function out = InitializeDefaultPythonEnv()
            cmd = 'conda create -y -n MRIToolkit python=3.10';
            out = ShellAndPython.ExecuteSH(cmd);
            if(contains(out,'conda deactivate'))
                disp('Success! Default python environment succesfully created.')
            else
                disp('Error creating default python environment.')
                disp(out)
            end
        end

    end
end