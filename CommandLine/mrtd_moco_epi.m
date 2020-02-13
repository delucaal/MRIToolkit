function mrtd_moco_epi(varargin)
    disp('mrtd_moco_epi');
    coptions = varargin;
    global MRIToolkit;
%     disp(varargin)

    file_in = GiveValueForName(coptions,'-help');
    if(~isempty(file_in) || isempty(varargin))
        help = 'This tools perform motion correction / Eddy currents correction / EPI distorion correction (optional) of diffusion MRI data';
        help = [help newline 'The gradient b-matrix is automatically reoriented.'];
        help = [help newline];
        help = [help newline 'usage: mrtd_moco_epi -nii file.nii -bval file.bval -bvec file.bvec -out corrected_file (other_options)'];
        help = [help newline];
        help = [help newline '-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x]' ....
                   ' 4=[x z y] =[y z x] 6=[z x y]'];
        help = [help newline '-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]'];
        help = [help newline '-epi: .nii file to perform EPI correction (T1 or T2 image)'];
        help = [help newline '-epi_constraint: allow deformation on specific axis. Input between quotas "1 1 1"'];
        help = [help newline '-epi_reg_mode: image to use for registration to structural. One between "fa" (default), "b0", "dwis"'];
        help = [help newline '-epi_normcorr: 0 or 1. Use normalized correlation in place of mutual information to drive the registration.'];        
        help = [help newline];
        fprintf(help);
        
        return
    end
    
    if(~isfield(MRIToolkit,'Elastix'))
        if(ismac)
           MRIToolkit.Elastix.ElastixCMD = fullfile(pwd,'MD_cor_E','macOSX64','elastix_Mac_64');
           MRIToolkit.Elastix.TransformixCMD = fullfile(pwd,'MD_cor_E','macOSX64','transformix_Mac_64');
           MRIToolkit.Elastix.Location = fullfile(pwd,'MD_cor_E','macOSX64');
%            disp('Locations:');
%            disp(MRIToolkit.Elastix.Location);
%            disp(MRIToolkit.Elastix.ElastixCMD);
%            disp('Current path');
%            disp(pwd);
%            disp('Code path');
%            disp(which('mrtd_moco_epi'));
        else
            disp('Not a mac');
        end
    else
        disp('Already initialized');
    end
    
    file_in = GiveValueForName(coptions,'-nii');
    if(isempty(file_in))
        error('Need to specify the target .nii file');
    end
    bval_file = GiveValueForName(coptions,'-bval');
    if(isempty(bval_file))
        error('Need to specify the target .bval file');
    end
    bvec_file = GiveValueForName(coptions,'-bvec');
    if(isempty(bvec_file))
        error('Need to specify the target .bvec file');
    end
    output_file = GiveValueForName(coptions,'-out');
    if(isempty(output_file))
        error('Need to specify the output name!');
    end
    grad_perm = GiveValueForName(coptions,'-grad_perm');
    if(isempty(grad_perm))
        grad_perm = [];
    else
        grad_perm = str2double(grad_perm);
    end
    grad_flip = GiveValueForName(coptions,'-grad_flip');
    if(isempty(grad_flip))
        grad_flip = [];
    else
        grad_flip = str2double(grad_flip);
    end
    epi_tgt = GiveValueForName(coptions,'-epi');
    if(isempty(epi_tgt))
        epi_tgt = [];
    end   
    epi_constraint = GiveValueForName(coptions,'-epi_constraint');
    if(isempty(epi_constraint))
        epi_constraint = [1 1 1];
    else
        pieces = strsplit(epi_constraint,' ');
        epi_constraint = [str2double(pieces{1}) str2double(pieces{2}) str2double(pieces{3})];
    end      
    epi_reg_mode = GiveValueForName(coptions,'-epi_reg_mode');
    if(isempty(epi_reg_mode))
        epi_reg_mode = 'fa';
    end      
    epi_normcorr = GiveValueForName(coptions,'-epi_normcorr');
    if(isempty(epi_normcorr))
        epi_normcorr = 1;
    else
        epi_normcorr = str2double(epi_normcorr);
    end      
    
    temp_mat_file = MRTD.QuickNiiBvalBvecToMat('nii_file',file_in,...
        'bval_file',bval_file,'bvec_file',bvec_file,'grad_perm',grad_perm,...
        'grad_flip',grad_flip);
    
    if(~isempty(epi_tgt))
        EDTI.PerformMocoEPI('mat_file',temp_mat_file,'epi_tgt',epi_tgt,'constraint_epi',epi_constraint,'epi_reg_mode',epi_reg_mode,'use_normcorr',epi_normcorr);
        mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',[temp_mat_file(1:end-4) '_MD_C_trafo.mat']);
    else
        EDTI.PerformMocoEPI('mat_file',temp_mat_file);
        mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',[temp_mat_file(1:end-4) '_MD_C_native.mat']);
    end
    
    [fp,fn,~] = fileparts(output_file);
    output_file = fullfile(fp,fn);
    EDTI.WriteNifti(mrt_data,[output_file '.nii']);
    
    fout = fopen([output_file '.bval'],'wt');
    for ij=1:length(mrt_data.bvals)
       fprintf(fout,'%f ',mrt_data.bvals(ij)); 
    end
    fclose(fout);
    
    fout = fopen([output_file '.bvec'],'wt');
    for ij=1:size(mrt_data.bvecs,2)
        for ik=1:size(mrt_data.bvecs,1)
           fprintf(fout,'%f ',mrt_data.bvecs(ik,ij)); 
        end
        fprintf(fout,'%s',newline);
    end
    fclose(fout);
    
    h = figure('Visible','off');

    filename = [output_file '.gif'];
    delay = 0.5;
    for n = 1:size(mrt_data.img,4)
      V = rot90(squeeze(mrt_data.img(round(end/2),:,:,n)));
      V2 = squeeze(mrt_data.img(:,:,round(end/2),:));
      subplot(1,2,1)
      imagesc(V,[0 prctile(V(:),95)]);
      axis image; % this ensures that getframe() returns a consistent size
      axis off;
      drawnow 
      subplot(1,2,2)
      imagesc(V2,[0 prctile(V2(:),95)]);
      axis image; % this ensures that getframe() returns a consistent size
      axis off;
      drawnow 
      colormap gray;
            
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',delay); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',delay); 
      end 
    end
    
    close all;
    
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
