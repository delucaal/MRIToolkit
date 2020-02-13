function mrtd_coordsys_fix(varargin)
    disp('mrtd_coordsys_fix');
    coptions = varargin;
%     disp(varargin)

    file_in = GiveValueForName(coptions,'-help');
    if(~isempty(file_in) || isempty(varargin))
        help = 'This tools can be used to check whether spatial dimension flip / permute are needed to properly process the given diffusion MRI data';
        help = [help newline];
        help = [help newline 'usage: mrtd_coordsys_fix -nii file.nii -bval file.bval -bvec file.bvec -out screenshot_file_no_extension (other_options)'];
        help = [help newline];
        help = [help newline '-perm: how to permute the spatial dimensions,e.g. "1 2 3" or "2 1 3" etc.'];
        help = [help newline '-flip: how to flip the sign of the diffusion gradients, e.g. "0 0 0" or "0 1 0" etc.'];
        help = [help newline];
        fprintf(help);
        
        return
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
    
    try
        orig_file = load_nii(file_in);
    catch
        reslice_nii(file_in,fullfile(tempdir,'reslice_nii.nii'));
        file_in = fullfile(tempdir,'reslice_nii.nii');
        orig_file = load_nii(file_in);
    end
        
    if(ndims(orig_file.img) == 4)
        orig_file.img = permute(orig_file.img,[2 1 3 4]);
    elseif(ndims(orig_file.img) == 3)
        orig_file.img = permute(orig_file.img,[2 1 3]);
    else
        error('Unexpected spatial dimensions');
    end
    orig_file.img = flip(orig_file.img,1);
    orig_file.img = flip(orig_file.img,2);
    
    EDTI.FlipPermuteSpatialDimensions('nii_file',file_in,'flip',spat_flip,...
        'permute',spat_perm,'output',fullfile(tempdir,'reslice_nii_test.nii'));
    
    data = EDTI.LoadNifti(fullfile(tempdir,'reslice_nii_test.nii'));
    
    h = figure('Visible','off','Color','k','InvertHardCopy','off');

    filename = [output_file '.png'];
    
    V2 = data.img(:,:,round(end/2),1);
    V1 = orig_file.img(:,:,round(end/2),1);
    
    subplot(2,2,1);
    imagesc(V1,[0 prctile(V1(:),95)]);
    axis image;
    axis off
    colormap gray;
    subplot(2,2,2);
    imagesc(V2,[0 prctile(V1(:),95)]);
    axis image;
    axis off
    colormap gray;

    V2 = rot90(squeeze(data.img(round(end/2),:,:,1)));
    V1 = rot90(squeeze(orig_file.img(round(end/2),:,:,1)));
    
    subplot(2,2,3);
    imagesc(V1,[0 prctile(V1(:),95)]);
    axis image;
    axis off
    colormap gray;
    subplot(2,2,4);
    imagesc(V2,[0 prctile(V1(:),95)]);
    axis image;
    axis off
    colormap gray;
    
    print(h,filename,'-dpng','-r300');
%     % Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     imwrite(imind,cm,filename); 
% 
%     close all;
    
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
