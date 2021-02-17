%% Compilation without CAT12
global MRIToolkit;
if(ispc)
    output_folder= getenv('USERPROFILE');
else
    output_folder= getenv('HOME');
end

mcc('-m','mrtcmd.m','-d',output_folder,'-v','-a',MRIToolkit.Elastix.Location,...
    '-a',MRIToolkit.dcm2niix,'-a',which('Grad_dirs_1024.txt'))


%% SPM based compilation
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_make_standalone.m 7483 2018-11-12 13:19:31Z guillaume $

global MRIToolkit;

%-Check startup.m
%--------------------------------------------------------------------------
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end

%-Input arguments
%--------------------------------------------------------------------------
if(exist('outdir','var') < 1)
%     if(ispc)
%         outdir = [getenv('USERPROFILE') filesep];
%     else
%         outdir = [getenv('HOME') filesep];
%     end
    outdir = uigetdir('Select the output folder');
end
contentsver = ''; 

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('Dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:numel(d)
    fi = spm_select('List',d{i},'.*tbx_cfg_.*\.m$');
    if ~isempty(fi)
        for ik=1:size(fi,1)
            if(strcmp(fi(ik,1),'.') < 1)
                ft = [ft(:); cellstr(fi(ik,:))];
            end
        end
    end
end
%-Create code to insert toolbox config
if isempty(ft)
    ftstr = '';
else
    ft = spm_file(ft,'basename');
    ftstr = sprintf('%s ', ft{:});
end
fprintf(fid,'values = {%s};\n', ftstr);
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
               fullfile(spm('Dir'),'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end
if ~isempty(contentsver)
    % Format: 'xxxx (SPMx) dd-mmm-yyyy'
    f = fileread(fullfile(spm('Dir'),'Contents.txt'));
    f = regexprep(f,'% Version \S+ \S+ \S+',['% Version ' contentsver]);
    fid = fopen(fullfile(spm('Dir'),'Contents.txt'),'w');
    fprintf(fid,'%s',f);
    fclose(fid);
end

%==========================================================================
%-Trim FieldTrip
%==========================================================================
d = fullfile(spm('Dir'),'external','fieldtrip','compat');
d = cellstr(spm_select('FPList',d,'dir'));
for i=1:numel(d)
    f = spm_file(d{i},'basename');
    nrmv = strncmp(f,'matlablt',8);
    if nrmv
        [dummy,I] = sort({f(9:end),version('-release')});
        nrmv = I(1) == 2;
    end
    if ~nrmv
        [sts, msg] = rmdir(d{i},'s');
    end
end

%==========================================================================
%-Compilation
%==========================================================================
Nopts = {'-p',fullfile(matlabroot,'toolbox','signal'),...
    '-p',fullfile(matlabroot,'toolbox','curvefit'),...
    '-p',fullfile(matlabroot,'toolbox','images'),...
    '-p',fullfile(matlabroot,'toolbox','optim'),...,...
    '-p',fullfile(matlabroot,'toolbox','stats')};
% if ~exist(Nopts{2},'dir'), Nopts = {}; end
Ropts = {};%{'-R','-singleCompThread'} ;
% if spm_check_version('matlab','8.4') >= 0
%     Ropts = [Ropts, {'-R','-softwareopengl'}];
% end
mcc('-m', 'mrtcmd.m', '-C', '-v',...
    '-o','mrtcmd',...
    '-d',outdir,...
    '-N',Nopts{:},...
    Ropts{:},...
    '-a',spm('Dir'),...
    '-a',MRIToolkit.Elastix.Location,...
    '-a',MRIToolkit.dcm2niix,...
    '-a',which('Grad_dirs_1024.txt') ...
    );


