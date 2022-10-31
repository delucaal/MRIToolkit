function mrtd_dti_dki(varargin)
disp('mrtd_dti_dki');
coptions = varargin;
if(isdeployed && length(varargin{1}) > 1)
    coptions = varargin{1};
end
%     disp(varargin)

if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
    help = 'This tool performs DTI or DKI fit of diffusion MRI data';
    help = [help newline];
    help = [help newline 'usage: mrtd_dti_dki -dki 0-1 -nii file.nii -bval file.bval -bvec file.bvec -out output_prefix (other_options)'];
    help = [help newline];
    help = [help newline '-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x]' ....
        ' 4=[x z y] =[y z x] 6=[z x y]'];
    help = [help newline '-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]'];
    help = [help newline '-estimator: one in "ols" (ordinary least squares),"wls" (weighted ls), "rekindle" (robust)'];
    help = [help newline '-constraint_:dki: stabilize DKI with the a DKI constraint (>0) - exclusive with the following'];
    help = [help newline '-mkcurve: stabilize DKI with the MK-Curve method (Fan et al.)'];
    help = [help newline];
    fprintf(help);
    
    return
end

nii_file = '';
bval_file = '';
bvec_file = '';
outfile = '';
grad_perm = [];
grad_flip = [];
estimator = 'ols';
mkcurve = 0;
dki = 0;
dki_constraints = [];

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
    elseif(strcmp(value,'-grad_perm'))
        grad_perm = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-grad_flip'))
        grad_flip = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-estimator'))
        estimator = coptions{input_id+1};
    elseif(strcmp(value,'-dki'))
        dki = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-mkcurve'))
        mkcurve = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-dki_constraints'))
        dki_constraints = str2double(coptions{input_id+1});
    end
    
end

if(isempty(nii_file))
    error('Missing mandatory argument -nii');
elseif(isempty(bval_file))
    error('Missing mandatory argument -bval');
elseif(isempty(bvec_file))
    error('Missing mandatory argument -bvec');
elseif(isempty(outfile))
    error('Missing mandatory argument -out');
end

while(true)
    dest_basename = fullfile(tempdir,['mrtd_' num2str(randi(500000))]);
    if(isempty(dir([dest_basename '*'])))
        break
    end
end

MRTQuant.b_Matrix_from_bval_bvec('bval_file',bval_file,'bvec_file',bvec_file,...
    'output',[dest_basename '.txt']);

while(true)
    temp_file = fullfile(tempdir,['tempnii_' num2str(randi(100000)) '.nii']);
    if(exist(temp_file,'file') < 1)
        break
    end
end

MRTQuant.ApplyRescaleSlope('nii_file',nii_file,'output',temp_file);
MRTQuant.ConformSpatialDimensions('nii_file',temp_file,'output',temp_file);

MRTQuant.PerformDTI_DKIFit('nii_file',temp_file,'txt_file',[dest_basename '.txt'],...
    'grad_perm',grad_perm,'grad_flip',grad_flip,'dki',dki,...
    'dki_constraints',dki_constraints,'output',[outfile '.mat'],...
    'mk_curve',mkcurve,'fit_mode',estimator);

MRTQuant.MatMetrics2Nii([outfile '.mat'],dki);

end
