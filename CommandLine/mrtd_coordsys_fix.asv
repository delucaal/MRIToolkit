function mrtd_coordsys_fix(varargin)
    disp('mrtd_coordsys_fix');
    coptions = varargin;
    % if(length(varargin{1}) > 1)
    %     coptions = varargin{1};
    % end
    %     disp(varargin)

    if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
        help = 'This tools can be used to check whether spatial dimension flip / permute are needed to properly process the given diffusion MRI data';
        help = [help newline];
        help = [help newline 'usage: mrtd_coordsys_fix -nii file.nii -out permuted_file.nii'];
        help = [help newline];
        help = [help newline '-perm: how to permute the spatial dimensions,e.g. "1 2 3" or "2 1 3" etc.'];
        help = [help newline '-flip: how to flip the spatial dimensions, e.g. "0 0 0" or "0 1 0" etc.'];
        help = [help newline '-auto: attempt auto detection (default 0) - overrides perm and flip'];
        help = [help newline];
        fprintf(help);
        
        return
    end
        
    file_in = GiveValueForName(coptions,'-nii');
    if(isempty(file_in))
        error('Need to specify the target .nii file');
    end
    auto_fix = GiveValueForName(coptions,'-auto');
    if(isempty(auto_fix))
        auto_fix = 0;
    end
    output_file = GiveValueForName(coptions,'-out');
    if(isempty(output_file))
        error('Need to specify the output name!');
    end
    if(auto_fix == 0)
        spat_perm = GiveValueForName(coptions,'-perm');
        if(isempty(spat_perm))
            spat_perm = [1 2 3];
        else
            C = strsplit(spat_perm,' ');
            spat_perm = [str2double(C(1)) str2double(C(2)) str2double(C(3))];
        end
        spat_flip = GiveValueForName(coptions,'-flip');
        if(isempty(spat_flip))
            spat_flip = [0 0 0];
        else
            C = strsplit(spat_flip,' ');
            spat_flip = [str2double(C(1)) str2double(C(2)) str2double(C(3))];
        end  

        MRTQuant.FlipPermuteSpatialDimensions('nii_file',file_in,'flip',spat_flip,...
            'permute',spat_perm,'output',output_file);
    else
        MRTQuant.ConformSpatialDimensions('nii_file',file_in,'output',output_file);
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
