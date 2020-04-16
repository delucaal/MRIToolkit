function mrtcmd(varargin)
    disp('mrtcmd');
    if(isempty(varargin) || isempty(varargin{1}) || strcmpi(varargin{1},'-help'))
        help = 'This is the command line tool of MRIToolkit!';
        help = [help newline];
        help = [help newline 'Available commands:'];
        help = [help newline 'mrtd_coordsys_fix - flip/permute spatial dimensions'];
        help = [help newline 'mrtd_test_gradients - detect flip/permute of the gradient directions'];
        help = [help newline 'mrtd_moco_epi - perform motion/EPI correction'];
        help = [help newline 'mrtd_dti_dki - DTI/DKI fit'];
        help = [help newline 'mrtd_deconv_fod - Perform spherical deconvolution'];
        help = [help newline];
        help = [help newline 'For help on usage, type mrtcmd command'];
        fprintf(help);
        
        disp('');
        which('mrtcmd')
        disp('');
        disp(ctfroot)
        return
    end
    warning('off');
    global MRIToolkit;
    MRIToolkit.Elastix.Location = fullfile(ctfroot,'/Users/albdl/Documents/ExploreDTI_Matlab/Source/MD_cor_E/macOSX64'); % Adjust per OS
    elastix_cmd = (fullfile(MRIToolkit.Elastix.Location,'elastix_Mac_64'));
    MRIToolkit.Elastix.ElastixCMD = elastix_cmd;
    transformix_cmd = (fullfile(MRIToolkit.Elastix.Location,'transformix_Mac_64'));
    MRIToolkit.Elastix.TransformixCMD = transformix_cmd;
    
    disp(MRIToolkit.Elastix);
    
    if(nargin > 1)
        options = varargin(2:end);
    else
        options = [];
    end
    cmd = varargin{1};
    if(strcmpi(cmd,'mrtd_coordsys_fix'))
        mrtd_coordsys_fix(options);
    elseif(strcmpi(cmd,'mrtd_test_gradients'))
        mrtd_test_gradients(options);
    elseif(strcmpi(cmd,'mrtd_moco_epi'))
        mrtd_moco_epi(options);
    elseif(strcmpi(cmd,'mrtd_dti_dki'))
        mrtd_dti_dki(options);
    elseif(strcmpi(cmd,'mrtd_deconv_fod'))
        mrtd_deconv_fod(options);
    else
        disp('Unknown command');
    end        
end
