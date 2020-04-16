function mrtd_test_gradients(varargin)
    coptions = varargin;
    if(length(varargin{1}) > 1)
        coptions = varargin{1};
    end
    %     disp(varargin)

    if(isempty(coptions) || isempty(coptions{1}) || strcmpi(coptions{1},'-help'))
        help = 'This tools can be used to check whether gradients flip / permute are needed to properly process the given diffusion MRI data';
        help = [help newline];
        help = [help newline 'usage: mrtd_test_gradients -nii file.nii -bval file.bval -bvec file.bvec -out screenshot_file_no_extension (other_options)'];
        help = [help newline];
        help = [help newline '-grad_perm: how to permute the diffusion gradients 1=[x y z] 2=[y x z] 3=[z y x]' ....
                   ' 4=[x z y] =[y z x] 6=[z x y]'];
        help = [help newline '-grad_flip: how to flip the sign of the diffusion gradients 1=[x y z] 2=[-x y z] 3=[x -y z] 4=[x y -z]'];
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
    
    temp_mat_file = MRTD.QuickNiiBvalBvecToMat('nii_file',file_in,...
        'bval_file',bval_file,'bvec_file',bvec_file,'grad_perm',grad_perm,...
        'grad_flip',grad_flip);
    FiberOrientation = load(temp_mat_file,'FE','FA');
    FiberOrientation.FA = FiberOrientation.FA/sqrt(3);
    
    h = figure('Visible','off','Color','k','InvertHardCopy','off');

    filename = [output_file '.png'];
    V = squeeze(FiberOrientation.FE(:,:,round(end/2),[2 1 3]));
    V2 = squeeze(FiberOrientation.FA(:,:,round(end/2)));
    
    subplot(1,2,1)
    imshow(abs(V2),[0 0.7]);
    axis image; % this ensures that getframe() returns a consistent size
    axis off;
    hold on
    for x=1:size(V,1)
        for y=1:size(V,2)
            Vec = squeeze(V(x,y,:));
            if(sum(Vec) == 0 || any(isnan(Vec)))
                continue
            end
            plot([y y+Vec(1)],[x x+Vec(2)],'-','Color',abs(Vec));
        end
    end
    drawnow
    
    V = squeeze(FiberOrientation.FE(round(end/2),:,:,[2 3 1]));
    V = permute(V,[2 1 3]);
    V = flip(V,1);
    V2 = squeeze(FiberOrientation.FA(round(end/2),:,:));
    V2 = permute(V2,[2 1]);
    V2 = flip(V2,1);
    subplot(1,2,2)
    imshow(abs(V2),[0 0.7]);
    axis image; % this ensures that getframe() returns a consistent size
    axis off;
    hold on
    for x=1:size(V,1)
        for y=1:size(V,2)
            Vec = squeeze(V(x,y,:));
            if(sum(Vec) == 0 || any(isnan(Vec)))
                continue
            end
            plot([y y+Vec(1)],[x x-Vec(2)],'-','Color',abs(Vec([1 3 2])));
        end
    end
    drawnow
    
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
