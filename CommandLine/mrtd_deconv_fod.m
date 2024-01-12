function mrtd_deconv_fod(varargin)
    disp('mrtd_deconv_fod');
    coptions = varargin;
    % if(length(varargin{1}) > 1)
    %     coptions = varargin{1};
    % end
    %     disp(varargin)

    if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
        help = 'This tool performs spherical deconvolution of (multi-shell) diffusion MRI data';
        help = [help newline];
        help = [help newline 'usage: mrtd_deconv_fod -method chosen_method -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file.nii (other_options)'];
        help = [help newline];
        help = [help newline '-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x]' ....
            ' 4=[x z y] =[y z x] 6=[z x y]'];
        help = [help newline '-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]'];
        help = [help newline '-method: one in "csd","mscsd","grl","grl_auto","mfod"'];
        help = [help newline '-t1_seg: T1 segmentation file derived from "fsl_anat" (mscsd only)'];
        help = [help newline '-lmax: spherical harmonics order'];
        help = [help newline '-mat: provide an ExploreDTI-compatible .MAT file. Overrides -nii -bval -bvec'];
        help = [help newline];
        help = [help newline 'GRL - mFOD specific options'];
        help = [help newline];
        help = [help newline '-aniso_rf_dti: Add an anisotropic RF using the DTI model. Specify the eigenvalues in ms2/um2 as "1.7 0.2 0.2"'];
        help = [help newline '-aniso_rf_dki: Add an anisotropic RF using the DKI model. Specify the eigenvalues in ms2/um2 and the mean kurtosis as "1.7 0.2 0.2 1" or "auto" to estimate from the data'];
        help = [help newline '-aniso_rf_noddi: Add an anisotropic RF using the NODDI model. Specify 6 parameters  "intra-cellular-volume free-diffusivity*10^9 watson-concentration isotropic-fraction isotropic-diffusivity b0-amplitude"'];
        help = [help newline '-iso_rf: Add an isotropic RF using the ADC model. Specify the diffusivity in um2/ms, as 0.7 for GM and 3 for CSF'];
        help = [help newline '-shell_weight: Inner-shell weighting factor (GRL only)'];
        help = [help newline];
        help = [help newline 'GRL default usage: mrtd_deconv_fod -method grl_auto -nii file.nii -bval file.bval -bvec file.bvec -out GRL_processing'];
        help = [help newline 'mFOD default usage: mrtd_deconv_fod -method mfod -nii file.nii -bval file.bval -bvec file.bvec -out mFOD_processing -aniso_rf_dki "auto" -aniso_rf_noddi "0.4 1.7 1 0 3 1" -iso_rf 3'];
        help = [help newline];
        fprintf(help);

        return
    end

nii_file = '';
bval_file = '';
bvec_file = '';
mat_file = '';
outfile = '';
grad_perm = 1;
grad_flip = 1;
deconv_method = 'csd';
t1_seg = '';
aniso_rf = {};
iso_rf = {};
lmax = 8;

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
    elseif(strcmp(value,'-method'))
        deconv_method = coptions{input_id+1};
    elseif(strcmp(value,'-lmax'))
        lmax = str2double(coptions{input_id+1});
    elseif(strcmp(value,'-t1_seg'))
        t1_seg = coptions{input_id+1};
    elseif(strcmp(value,'-aniso_rf_dti'))
        aniso_rf(end+1) = {{'dti',coptions{input_id+1}}};
    elseif(strcmp(value,'-aniso_rf_dki'))
        aniso_rf(end+1) = {{'dki',coptions{input_id+1}}};
    elseif(strcmp(value,'-aniso_rf_noddi'))
        aniso_rf(end+1) = {{'aniso_rf_noddi',coptions{input_id+1}}};
    elseif(strcmp(value,'-iso_rf'))
        iso_rf(end+1) = {coptions{input_id+1}};
    elseif(strcmp(value,'-mat'))
        mat_file = coptions{input_id+1};
    end
    
end

if(isempty(deconv_method))
    error('Missing mandatory argument -method');
end

if(isempty(mat_file))
    if(isempty(nii_file))
        error('Missing mandatory argument -nii');
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

if(strcmp(deconv_method,'csd') || strcmp(deconv_method,'mscsd'))
    if(strcmp(deconv_method,'mscsd'))
        mscsd = 1;
    else
        mscsd = 0;
    end
    EDTI.PerformCSD('mat_file',temp_mat_file,'mscsd',mscsd,...
        'T1seg',t1_seg,'output',outfile,'lmax',lmax);

elseif(strcmp(deconv_method,'grl') || strcmp(deconv_method,'mfod'))
    % Converts a .MAT to the MRIToolkit format
    mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',temp_mat_file);
    other_props = load(temp_mat_file,'FA');
    % Prepare the spherical deconvolution class on this dataset
    SD = MRTTrack('data',mrt_data);
    % Add an anisotropic response function using the DKI model
    for rf_id=1:length(aniso_rf)
        requested_rf = aniso_rf{rf_id};
        eigvals = requested_rf{2};
        if(contains(requested_rf{1},'dki') > 0 && contains(eigvals,'auto'))
            [eigvals,isok] = MRTTrack.Eigenval_IsoK_WM_FromData(mrt_data,other_props.FA,~isnan(other_props.FA));
            SD.AddAnisotropicRF_DKI(eigvals,isok); % WM
        else
            J = strsplit(eigvals,' ');
            eigvals = [];
            for ikk=1:length(J)
               eigvals = [eigvals str2double(J{ikk})]; 
            end
            if(contains(requested_rf{1},'dti'))
                isok = 0;
                eigvals = eigvals(1:3)*1e-3;
                SD.AddAnisotropicRF_DKI(eigvals,isok); % WM
            elseif(contains(requested_rf{1},'dki'))
                isok = eigvals(4);
                eigvals = eigvals(1:3)*1e-3;
                SD.AddAnisotropicRF_DKI(eigvals,isok); % WM
            elseif(contains(requested_rf{1},'noddi'))
                SD.AddAnisotropicRF_NODDI(eigvals);
            end
        end
    end
    
    for rf_id=1:length(iso_rf)
        SD.AddIsotropicRF(str2double(iso_rf{rf_id})*1e-3);  
    end
    
    if(strcmp(deconv_method,'grl'))
        SD.setInnerShellWeighting(0.2); % Do not loose angular resolution due to the lower shells
        SD.AutomaticDRLDamping(); % see Dell'Acqua 2013
        SD.setDeconvMethod('dRL'); % damped Richardson Lucy
    elseif(strcmp(deconv_method,'mfod'))
        SD.setInnerShellWeighting(1.0); % in mFOD this parameter has not been investigated
        SD.setDeconvMethod('L2LSQ');
    end
    
    GRL_Results = SD.PerformDeconv();
    MRTTrack.SaveOutputToNii(SD,GRL_Results,outfile);
elseif(strcmp(deconv_method,'grl_auto'))
    MRTTrack.PerformGRL('mat_file',temp_mat_file,'output',outfile); 
else
    error(['Unknown deconvolution method ' method]);
end

end

