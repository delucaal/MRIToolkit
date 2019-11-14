%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function MRIT_DTIFit(varargin)
    fit_dki = 0;
    fit_mode = 0;
    constraint = 0;
    perm_mode = 1;
    flip_mode = 1;
        
    if(nargin == 0)
        disp('Usage: MRIT_DTIFit target(nii or folder) -dki -ols[or -wls -robust] -pos_dki -grad_perm 1 -grad_flip 1');
        return
    end
    
    target = varargin{1};
    
    for arg_id=2:length(varargin)
        if(strcmpi(varargin{arg_id},'dki'))
            fit_dki = 1;
        elseif(strcmpi(varargin{arg_id},'pos_dki'))
            constraint = 1;
        elseif(strcmpi(varargin{arg_id},'ols'))
            fit_mode = 0;
        elseif(strcmpi(varargin{arg_id},'wls'))
            fit_mode = 1;
        elseif(strcmpi(varargin{arg_id},'robust'))
            fit_mode = 2;
        elseif(strcmpi(varargin{arg_id},'grad_perm'))
            perm_mode = varargin{arg_id};
        elseif(strcmpi(varargin{arg_id},'grad_flip'))
            flip_mode = varargin{arg_id};
        end
    end
    
    if(fit_mode == 0)
        fit_mode = 'ols';
    elseif(fit_mode == 1)
        fit_mode = 'wls';
    elseif(fit_mode == 2)
        fit_mode = 'rekindle';
    end
    
    if(constraint == 1)
        constraint = [0 Inf];
    else
        constraint = [-Inf Inf];
    end
    
    if(isfolder(target))
        nii_files = dir(fullfile(target,'*.nii'));
        for file_id=1:length(nii_files)
            try
                ExploreDTI_Interface.PerformDTI_DKIFit('nii_file',fullfile(target,nii_files(file_id).name),'grad_perm',perm_mode,'grad_flip',flip_mode,'fit_mode',fit_mode,'dki_constraints',constraint,'dki',fit_dki);
            catch
                disp(['Error trying to fit ' nii_files(file_id).name]);
            end
        end
    else        
        ExploreDTI_Interface.PerformDTI_DKIFit('nii_file',target,'grad_perm',perm_mode,'grad_flip',flip_mode,'fit_mode',fit_mode,'dki_constraints',constraint,'dki',fit_dki);
    end

end