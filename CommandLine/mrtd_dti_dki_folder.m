function mrtd_dti_dki_folder(varargin)
disp('mrtd_dti_dki_folder');
coptions = varargin;
if(length(varargin{1}) > 1)
    coptions = varargin{1};
end
%     disp(varargin)

if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
    help = 'This tool performs DTI or DKI fit of diffusion MRI data in an entire folder';
    help = [help newline];
    help = [help newline 'usage: mrtd_dti_dki_folder -dki 0-1 -folder target_folder (other_options)'];
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

tgt_folder = '';
grad_perm = [];
grad_flip = [];
estimator = 'ols';
mkcurve = 0;
dki = 0;
dki_constraints = [];

for input_id =1:2:length(coptions)
    value = coptions{input_id};
    if(strcmp(value,'-folder'))
        tgt_folder = coptions{input_id+1};
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

if(isempty(tgt_folder))
    error('Missing mandatory argument -folder');
end

files2fit = dir(fullfile(tgt_folder,'*.bval'));
for fileid=1:length(files2fit)
    basename = fullfile(files2fit(fileid).folder,files2fit(fileid).name(1:end-5));
    try
        bval = [basename '.bval'];
        bvec = [basename '.bvec'];
        nii = [basename '.nii'];
        mrtd_dti_dki('-nii',nii,'-bval',bval,'-bvec',bvec,...
          '-out',basename,'-dki',dki,'-mkcurve',mkcurve,'-dki_constraints',...
          dki_constraints,'-estimator',estimator,'-grad_perm',grad_perm,...
          '-grad_flip',grad_flip);
    catch
        disp(['Error while fitting ' basename]);
    end
end

end
