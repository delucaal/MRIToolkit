function mrtd_track(varargin)
    disp('mrtd_track');
    coptions = varargin;
    % if(nargin > 0 && length(varargin{1}) > 1)
    %     coptions = varargin{1};
    % else
    %     coptions = {};
    % end
    %     disp(varargin)

    if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
        help = 'This tool performs (deterministic) fiber tractography based on either DTI or FOD.';
        help = [help newline];
        help = [help newline 'DTI-based usage: mrtd_track -nii file.nii -bval file.bval -bvec file.bvec -out track.mat (other_options)'];
        help = [help newline 'FOD-based usage: mrtd_track -nii file.nii -bval file.bval -bvec file.bvec  -fod file.nii -out track.mat (other_options)'];
        help = [help newline];
        help = [help newline '-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x]' ....
            ' 4=[x z y] =[y z x] 6=[z x y] (DTI only)'];
        help = [help newline '-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z] (DTI only)'];
        help = [help newline '-seed_res: the seeding resolution in mm (default 2)'];
        help = [help newline '-step: the step size in mm (default 1)'];
        help = [help newline '-angle: angle threshod in degrees (default 30)'];
        help = [help newline '-fod_thresh: amplitude threshold of the FOD peaks (default 0.1)'];
        help = [help newline '-seed_mask: a .nii file containing a seeding mask'];
        help = [help newline '-mat: provide an ExploreDTI-compatible .MAT file. Overrides -nii -bval -bvec'];
        help = [help newline];
        fprintf(help);

        return
    end

nii_file = '';
bval_file = '';
bvec_file = '';
outfile = '';
fod = '';
grad_perm = [];
grad_flip = [];
seed_res = [2 2 2];
step_size = 1;
angle_thresh = 30;
fod_thresh = 0.1;
mat_file = '';
seed_mask_file = '';

for input_id =1:2:length(coptions)
    value = coptions{input_id};
    if(strcmp(value,'-nii'))
        nii_file = coptions{input_id+1};
    elseif(strcmp(value,'-bval'))
        bval_file = coptions{input_id+1};
    elseif(strcmp(value,'-bvec'))
        bvec_file = coptions{input_id+1};
    elseif(strcmp(value,'-out'))
        outfile = coptions{input_id+1};
    elseif(strcmp(value,'-fod'))
        fod = (coptions{input_id+1});
    elseif(strcmp(value,'-seed_res'))
        seed_res = str2double(coptions{input_id+1});
        seed_res = [seed_res seed_res seed_res];
    elseif(strcmp(value,'-step'))
        step_size = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-angle'))
        angle_thresh = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-fod_thresh'))
        fod_thresh = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-grad_perm'))
        grad_perm = coptions{input_id+1};
    elseif(strcmp(value,'-grad_flip'))
        grad_flip = coptions{input_id+1};
    elseif(strcmp(value,'-mat'))
        mat_file = coptions{input_id+1};
    elseif(strcmp(value,'-seed_mask'))
        seed_mask_file = coptions{input_id+1};
    end
    
end

if(isempty(outfile))
    error('Missing mandatory argument -out');
end

if(isempty(mat_file))
    if(isempty(nii_file))
        error('Missing mandatory argument -nii or -fod');
    elseif(isempty(bval_file))
        error('Missing mandatory argument -bval');
    elseif(isempty(bvec_file))
        error('Missing mandatory argument -bvec');
    end
    temp_mat_file = MRTD.QuickNiiBvalBvecToMat('nii_file',nii_file,...
    'bval_file',bval_file,'bvec_file',bvec_file,'grad_perm',grad_perm,...
    'grad_flip',grad_flip);
else
    temp_mat_file = mat_file;
end

if(contains(outfile,'.mat') < 1)
    outfile = [outfile '.mat'];
end

if(isempty(fod))
    disp('Performing DTI based tractography');
    MRTTrack.PerformDTIBased_FiberTracking('mat_file',temp_mat_file,...
        'SeedPointRes', seed_res, 'AngleThresh', angle_thresh,...
        'StepSize', step_size, 'SeedMask', seed_mask_file,...
        'output',outfile);    
else
    disp('Performing FOD based tractography');
    MRTTrack.PerformFODBased_FiberTracking('mat_file',temp_mat_file,...
        'fod_file',fod,'SeedPointRes', seed_res, 'AngleThresh', angle_thresh,...
        'StepSize', step_size, 'SeedMask', seed_mask_file,...
        'FODThresh',fod_thresh,'output',outfile);
end

MRTTrack.ConvertMatTractography2TCK('mat_file',outfile,'output',strrep(outfile,'.mat','.tck'));

end
