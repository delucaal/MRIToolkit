
%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% A. De Luca
% This function adds the specified toolboxes to the path. If no input is
% specified all toolboxes will be added
% 29/03/2018: creation - v1.0
% 23/09/2018: addition of ThirdParty, OptimizationMethods, Relaxometry
% 09/11/2019: Change of name - added the - ExploreDTIInterface
% 01/02/2020: Added the MRT class and CommandLineTools
% 05/05/2020: Added suppot for automatic version control
% 13/08/2021: v1.5
% 12/01/2024: v1.6
% 18/12/2024: v1.7
function MRIToolkitInit(SelectedToolboxes)
    run_folder = mfilename('fullpath');
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
    addpath(fullfile(run_folder,'init'));

    global MRIToolkit;
            
    MRIToolkit.version = read_git_version();
    MRIToolkit.EnforceNiiGz = false;

    try
       MRIToolkitDefineLocalVars(); 
    catch
        disp('I cannot find a configuration file. Please define one');
    end
    
    available_toolboxes = {
        {'NiftiIO_basic','Manages basic Nifti input/output',true},...
        {'ImageRegistrations','Elastix based registration utils',true},...
        {'ThirdParty','Third party utilities',true},...
        {'ExploreDTIInterface','ExploreDTI powered library',true},...
        {'Core','MRIToolkit core implementation',true},...
        {'Neuroimaging','Common neuro pipelines',true},...
        {'ReportMaker','Utilities to generate HTML-based reports of entire studies',true},...
        {'CommandLine','Command line tools (Executables)',true},...
        {'Legacy','Backward compatibility (deprecated)',true},...
        }; % Folder, Description, Load default
    
    if(nargin > 0 && ~iscell(SelectedToolboxes))
        disp('Usage 1: MRIToolkitInit() -> load all toolboxes');
        disp('Usage 2: MRIToolkitInit({''Identifier1'',...)');
        disp('Possible identifiers:');
        for ij=1:length(available_toolboxes)
           disp([available_toolboxes{ij}{1} ': ' available_toolboxes{ij}{2}]);
        end
        return
    end
    
    run_folder = get_executed_file_path();
    MRIToolkit.RootFolder = run_folder;

    % Add all toolboxes if no input is specified
    if(nargin < 1)
        SelectedToolboxes = available_toolboxes;
        override_default = false;
    else
        UserSpecifiedToolboxes = SelectedToolboxes;
        override_default = true;
        warning('off');
        disp('Deactivating previously initialized toolboxes (if any)');
        for tool_id=1:length(available_toolboxes)
            run(fullfile(run_folder,available_toolboxes{tool_id}{1},'mk_deinit.m'));
        end
        warning('on');
        SelectedToolboxes = cell(length(UserSpecifiedToolboxes),1);
        for tool_id=1:length(SelectedToolboxes)
            toolbox_found = 0;
            for av_id=1:length(available_toolboxes)
               if(strcmpi(UserSpecifiedToolboxes{tool_id},available_toolboxes{av_id}{1}))
                  toolbox_found = av_id;
                  break
               end
            end
           if(toolbox_found > 0)
              SelectedToolboxes{tool_id} = available_toolboxes{toolbox_found};
           else
              disp(['The specified toolbox (' UserSpecifiedToolboxes{tool_id} ') does not exist. Quitting.']);
              return
           end
        end
    end
    
    for tool_id=1:length(SelectedToolboxes)
        % Call the init function of each toolbox
        if(SelectedToolboxes{tool_id}{3} || override_default)
            disp(['Adding ' SelectedToolboxes{tool_id}{1} ' (' SelectedToolboxes{tool_id}{2} ') to the path']);
            run(fullfile(run_folder,SelectedToolboxes{tool_id}{1},'mk_init.m'));
        end
    end
    
end

function git_version = read_git_version()

    mrt_folder = which('MRIToolkitInit.m');
    fp = fileparts(mrt_folder);
    git_v = fullfile(fp,'.git','refs','heads','master');
    if(exist(git_v,'file') > 0)
        git_commit = fopen(git_v,'rt');
        git_version = fgetl(git_commit);
        fclose(git_commit);
    elseif(isdeployed)
        git_version = 'v1.7-cmdline';
    else
        git_version = 'v1.7-nongit';
    end
end
