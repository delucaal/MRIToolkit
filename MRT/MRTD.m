% This class implements methods for the MRIToolkit diffusion pipeline
% Differently from other classes (as EDTI), the input / output of this
% class are Niftis and other exchangeble formats (in addition to
% ExploreDTI-compatible .mat files)
classdef MRTD < handle
    methods(Static)
        
        % Create a .MAT from .nii and .bval/.bvec for internal processing
        % in a temporary directory. Input arguments:
        % nii_file: the target .nii file
        % bval_file: the compation .bval fie
        % bvec_file: the companion .bvec file
        % grad_perm: optional, how to permute the diffusion gradients
        % grad_flip: optional, how to flip the sign of the gradients
        function temp_mat_location = QuickNiiBvalBvecToMat(varargin)
            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Need to specify the target file');
            end
            
            bval_file = GiveValueForName(coptions,'bval_file');
            if(isempty(bval_file))
                bval_file = [file_in(1:end-4) '.bval'];
            end
            
            bvec_file = GiveValueForName(coptions,'bvec_file');
            if(isempty(bvec_file))
                bvec_file = [file_in(1:end-4) '.bvec'];
            end
            
            option_value = GiveValueForName(coptions,'grad_perm');
            if(~isempty(option_value))
                perm = option_value;
            else
                perm = [];
            end
            option_value = GiveValueForName(coptions,'grad_flip');
            if(~isempty(option_value))
                flip = option_value;
            else
                flip = [];
            end
            
            while(true)
                dest_basename = fullfile(tempdir,['mrtd_' num2str(randi(500000))]);
                if(isempty(dir([dest_basename '*'])))
                    break
                end
            end
            
            EDTI.b_Matrix_from_bval_bvec('bval_file',bval_file,'bvec_file',bvec_file,...
                'output',[dest_basename '.txt']);
            
            if(~isempty(perm) && ~isempty(flip))
                EDTI.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
                    'output',[dest_basename '.mat'],'grad_perm',perm,'grad_flip',flip);
            else
                EDTI.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
                    'output',[dest_basename '.mat']);
            end
            
            temp_mat_location = [dest_basename '.mat'];
            
        end
        
    end
end

% Helper: finds a parameter by name when using varargin
function value = GiveValueForName(coptions,name)
value = [];
for ij=1:2:length(coptions)
    if(strcmpi(coptions{ij},name))
        value = coptions{ij+1};
        return
    end
end
end
