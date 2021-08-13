function pi1c(varargin)
    disp('MRIToolkit - processing in one command (pi1c)');
    if(isempty(varargin) || isempty(varargin{1}) || strcmpi(varargin{1},'-help'))
        help = 'This is the command line tool of MRIToolkit!';
        help = [help newline];
        help = [help newline 'Available commands:'];
        help = [help newline 'mrtd_coordsys_check - check whether a flip/permute of the spatial dimensions is needed'];
        help = [help newline 'mrtd_coordsys_fix - flip/permute the spatial dimensions as needed'];
        help = [help newline 'mrtd_test_gradients - detect flip/permute of the gradient directions'];
        help = [help newline 'mrtd_preproc - perform dMRI preprocessing'];
        help = [help newline 'mrtd_dti_dki - DTI/DKI fit'];
        help = [help newline 'mrtd_deconv_fod - Perform spherical deconvolution'];
        help = [help newline 'mrtd_track - Fiber tractography'];
        help = [help newline 'mrtd_neuro - Popular neuroimaging pipelines'];
        help = [help newline];
        help = [help newline 'For help on usage, type mrtcmd command'];
        fprintf(help);
                
        return
    end
    warning('off');
    global MRIToolkit;
    if(isdeployed)
        if(ismac)
            MRIToolkit.Elastix.Location = fullfile(ctfroot,'/Users/albdl/Documents/ExploreDTI_Matlab/Source/MD_cor_E/macOSX64'); % Adjust per OS
            elastix_cmd = (fullfile(MRIToolkit.Elastix.Location,'elastix_Mac_64'));
            transformix_cmd = (fullfile(MRIToolkit.Elastix.Location,'transformix_Mac_64'));
            MRIToolkit.dcm2niix = fullfile(ctfroot,'/Applications/dcm2niix-master/bin/dcm2niix');
            MRIToolkit.spm_path = fullfile(ctfroot,'/Volumes/KINGSTON/Work/spm12/');
            MRIToolkit.noddi_path = fullfile(ctfroot,'/Volumes/KINGSTON/Work/NODDI_toolbox_v1.01/');
        elseif(isunix)
            MRIToolkit.Elastix.Location = fullfile(ctfroot,'/Users/albdl/Documents/ExploreDTI_Matlab/Source/MD_cor_E/linux64'); % Adjust per OS
            elastix_cmd = (fullfile(MRIToolkit.Elastix.Location,'elastix_L_64'));
            transformix_cmd = (fullfile(MRIToolkit.Elastix.Location,'transformix_L_64'));
        elseif(ispc)
            MRIToolkit.Elastix.Location = fullfile(ctfroot,'/Users/albdl/Documents/ExploreDTI_Matlab/Source/MD_cor_E/win64'); % Adjust per OS
            elastix_cmd = (fullfile(MRIToolkit.Elastix.Location,'elastix_64'));
            transformix_cmd = (fullfile(MRIToolkit.Elastix.Location,'transformix_64'));
        end
        MRIToolkit.Elastix.ElastixCMD = elastix_cmd;
        MRIToolkit.Elastix.TransformixCMD = transformix_cmd;
    end
    MRIToolkit.version = 'mrtcmd-1.1';
        
    if(nargin > 1)
        options = varargin(2:end);
    else
        options = [];
    end
    cmd = varargin{1};
    if(strcmpi(cmd,'mrtd_coordsys_check'))
        mrtd_coordsys_check(options);
    elseif(strcmpi(cmd,'mrtd_coordsys_fix'))
        mrtd_coordsys_fix(options);
    elseif(strcmpi(cmd,'mrtd_test_gradients'))
        mrtd_test_gradients(options);
    elseif(strcmpi(cmd,'mrtd_preproc'))
        mrtd_preproc(options);
    elseif(strcmpi(cmd,'mrtd_dti_dki'))
        mrtd_dti_dki(options);
    elseif(strcmpi(cmd,'mrtd_deconv_fod'))
        mrtd_deconv_fod(options);
    elseif(strcmpi(cmd,'mrtd_track'))
        mrtd_track(options);
    elseif(strcmpi(cmd,'mrtd_neuro'))
        mrtd_neuro(options);
    else
        disp('Unknown command');
    end        
end
