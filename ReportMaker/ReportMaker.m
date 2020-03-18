
classdef ReportMaker < handle
   properties
       datafolder
       folders2scan
       images2scan
       output_folder
       current_html
   end
   
   methods
      
       function out = ReportMaker()
           out.folders2scan = {};
           out.images2scan = {};
           out.output_folder = '';
           out.current_html = '';
           out.datafolder = ''; 
       end
       
       function me = AddFolder(me,folder)
            me.folders2scan(end+1) = {folder};
       end
       
       function me = AddDataFolder(me,folder)
           me.datafolder = folder;
       end
       
       function me = AddFolders(me,root_folder)
          folders = dir(root_folder);
          for foid=1:length(folders)
             if(strcmpi(folders(foid).name(1),'.') > 0)
                 continue
             end
             me.AddFolder(fullfile(root_folder,folders(foid).name));
          end
       end
       
       function me = AddImage2Scan(me,label,rel_path_to_file,search_string,axial_coronal_sagittal,link_to_viewer,make_gif,multi_slice,colormap,overlay_on,stride)
            img.label = label;
            img.rel_path_to_file = rel_path_to_file;
            img.search_string = search_string;
            img.axial_coronal_sagittal = axial_coronal_sagittal;
            img.link_to_viewer = link_to_viewer;
            img.make_gif = make_gif;
            img.multi_slice = multi_slice;
            img.colormap = colormap;
            img.overlay_on = overlay_on;
            img.stride = stride;
            me.images2scan(end+1) = {img};
       end
       
       function me = Process(me,out_folder)
            me.output_folder = out_folder;
            if(exist(me.output_folder,'dir') < 1)
                mkdir(me.output_folder);
            end
            html_code = '<html>';
            html_code = concat_html(html_code,'<body bgcolor="000000">');
            html_code = concat_html(html_code,'<table style="color:white">');
            html_code = concat_html(html_code,'<tr>');
            html_code = concat_html(html_code,'<th>Dataset</th>');
            for ij=1:length(me.images2scan)
                html_code = concat_html(html_code,['<th>' me.images2scan{ij}.label '</th>']);
            end
            html_code = concat_html(html_code,'</tr>');

            f=figure('Visible','off');
            set(gcf,'color','k','inverthardcopy','off')

            cdir = pwd;
            cd(out_folder);
            
            if(exist('Images','dir') < 1)
                mkdir('Images');
            end

            if(exist('Niftis','dir') < 1)
                mkdir('Niftis');
            end
            
            width_in_px = 300;
            ax_or_coronal = 1;
            try
                for subj_id=1:length(me.folders2scan)
                    the_folder = me.folders2scan{subj_id};
                    spf = strsplit(the_folder,'/');
                    if(length(spf) == 1)
                        spf = strsplit(the_folder,'\');
                    end
                    spf = spf{end};
                    html_code = concat_html(html_code,'<tr>');
                    html_code = concat_html(html_code,['<th>' spf '</th>']);
                    for ij=1:length(me.images2scan)
                       cel = me.images2scan{ij};
                       cmap = cel.colormap;
                       if(ischar(cmap))
                           cmap = colormap(cmap);
                       end
                       
                       html_code = concat_html(html_code,'<th>');
                      
                       file = dir(fullfile(me.datafolder,the_folder,cel.rel_path_to_file,cel.search_string));
                       if(isempty(file) || length(file) > 1)
                           warning('CHECK 1');
                           html_code = concat_html(html_code,'</th>');
                           continue
                       end

                       [~,fn,~] = fileparts(file.name);
                       if(cel.make_gif == 0)
                           img_name = ['Images/' spf '_' fn '.png'];
                       else
                           img_name = ['Images/' spf '_' fn '.gif'];
                       end
                       nii_name = fullfile(file.folder,file.name);
                       if(cel.make_gif == 0)
                          if(cel.multi_slice == 0)
                              if(cel.link_to_viewer == 0)
                                  html_code = concat_html(html_code,html_img_tag(img_name,nii_name,ax_or_coronal,width_in_px,cmap));
                              else
                                  overlay_img = 0;
                                  if(cel.overlay_on ~= 0)
                                     tfile = dir(fullfile(me.datafolder,the_folder, me.images2scan{cel.overlay_on}.rel_path_to_file, me.images2scan{cel.overlay_on}.search_string));
                                     [~,tfn,~] = fileparts(tfile.name);
                                     overlay_img = ['Niftis/' spf '_' tfn '.nii'];
                                  end
                                  html_code = concat_html(html_code,html_img_tag_3dviewer(img_name,nii_name,ax_or_coronal,width_in_px,'Niftis',cmap,overlay_img));
                              end
                          else
                              Vol3D = EDTI.LoadNifti(nii_name);
                              Vol3D = Vol3D.img;
                              html_code = concat_html(html_code,html_img_tag_multislice(img_name,Vol3D,width_in_px,cmap));
                              clear Vol3D;
                          end
                       else
                          if(cel.link_to_viewer == 0)
                              html_code = concat_html(html_code,make_gif(img_name,nii_name,ax_or_coronal,width_in_px,cmap,cel.stride));
                          else
                              overlay_img = 0;
                              if(cel.overlay_on ~= 0)
                                 tfile = dir(fullfile(me.datafolder,the_folder, me.images2scan{cel.overlay_on}.rel_path_to_file, me.images2scan{cel.overlay_on}.search_string));
                                 [~,tfn,~] = fileparts(tfile.name);
                                 overlay_img = ['Niftis/' spf '_' tfn '.nii'];
                              end
                              html_code = concat_html(html_code,make_gif_3dviewer(img_name,nii_name,ax_or_coronal,width_in_px,'Niftis',cmap,overlay_img,cel.stride));
                          end
                       end

                       html_code = concat_html(html_code,'</th>');
                    end
                    html_code = concat_html(html_code,'</tr>');
                end

                cd(cdir);

                html_code = concat_html(html_code,'</table>');
                html_code = concat_html(html_code,'</body>');
                html_code = concat_html(html_code,'</html>');

                me.current_html = html_code;
                write_html(html_code,fullfile(out_folder,'Index.html'));
                copyfile(which('papaya.js'),out_folder);
                copyfile(which('papaya.html'),out_folder);
            catch err
                disp(err.message);
                cd(cdir);
            end
            close(f);
       end
       
   end
   
end

function out_string = concat_html(in_string,piece)
    out_string = [in_string newline piece];
end

function write_html(html_code,file)
    fout = fopen(file,'wt');
    pieces = strsplit(html_code,'\n');
    for ij=1:length(pieces)
        fprintf(fout,'%s\r\n',strrep(pieces{ij},'\n',''));
    end
    fclose(fout);
end

function html_code = html_img_tag_multislice(img_name,Vol3D,width_in_px)
        if(nargin < 2)
            res = '-r150';
        else
            res = ['-r' num2str(res)];
        end
        if(nargin < 3)
            width_in_px = 300;
        end

        html_code = ['<a href="' img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
        E_DTI_imagescn_one_slice(flip(permute(Vol3D,[3 1 2]),1),m,VD); colormap hot;
        print(gcf,img_name,'-dpng',res);
end

function html_code = html_img_tag(img_name,nii_name,ax_or_coronal,width_in_px,cmap)
    img = EDTI.LoadNifti(nii_name);
    img = img.img;
    img = single(img);
    img = img/prctile(img(:),99);
    if(nargin < 4)
        width_in_px = 300;
        cmap = colormap;
    end
    if(ax_or_coronal == 1)
        imagesc(squeeze(img(:,:,round(end/2),:)),[0 1]);axis off;axis image;
        colormap(cmap);
    end
    html_code = ['<a href="' img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    print(gcf,img_name,'-dpng','-r150');
end

function html_code = html_img_tag_3dviewer(img_name,nii_name,ax_or_coronal,width_in_px,nifti_folder,cmap,overlay_img)
    img = EDTI.LoadNifti(nii_name);
    img = img.img;
    img = single(img);
    img = img/prctile(img(:),99);
    if(nargin < 4)
        width_in_px = 300;
        cmap = colormap;
    end
    if(ax_or_coronal == 1)
        imagesc(squeeze(img(:,:,round(end/2),:)),[0 1]);axis off;axis image;
        colormap(cmap)
        overlay_img = 0;
    end

    [out_folder,out_name] = fileparts(img_name);
    new_nii_name = fullfile(nifti_folder,[out_name '.nii']);
    copyfile(nii_name,new_nii_name);
    if(overlay_img == 0)
        html_code = ['<a href="papaya.html?image_uri=' urlencode(strrep(new_nii_name,'\','/')) '" '... 
                img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    else
        html_code = ['<a href="papaya.html?overlay_uri=' urlencode(strrep(new_nii_name,'\','/')) ...
            '&image_uri=' urlencode(strrep(overlay_img,'\','/')) '" '...
                img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    end            
    print(gcf,img_name,'-dpng','-r150');
end

function html_code = make_gif(img_name,nii_file,ax_or_coronal,width_in_px,cmap,stride)
    if(nargin < 4)
        width_in_px = 300;
        cmap = colormap;
        stride = 1;
    end
    img = EDTI.LoadNifti(nii_file);
    img = img.img;
    img = single(img);
    img = img/prctile(img(:),99);
    if(ax_or_coronal == 1)
      % Capture the plot as an image 
      if(ndims(img) == 3)
          for z=1:stride:size(img,3)
              imagesc(img(:,:,z),[0 1]);
              axis image;
              axis off;
              colormap(cmap)
              frame = getframe(gcf); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256); 
              % Write to the GIF File 
              if z == 1 
                  imwrite(imind,cm,img_name,'gif', 'Loopcount',inf); 
              else 
                  imwrite(imind,cm,img_name,'gif','WriteMode','append','DelayTime',0.1); 
              end  
          end
      elseif(ndims(img) == 4)
          for z=1:stride:size(img,4)
              imagesc(img(:,:,round(end/2),z),[0 1]);
              colormap(cmap);
              axis image;
              axis off;
              frame = getframe(gcf); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256); 
              % Write to the GIF File 
              if z == 1 
                  imwrite(imind,cm,img_name,'gif', 'Loopcount',inf); 
              else 
                  imwrite(imind,cm,img_name,'gif','WriteMode','append','DelayTime',0.1); 
              end           
          end
      end
    end
    
    html_code = ['<a href="' img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    
end

function html_code = make_gif_3dviewer(img_name,nii_file,ax_or_coronal,width_in_px,nifti_folder,cmap,overlay_img,stride)
    if(nargin < 4)
        width_in_px = 300;
        cmap = colormap;
        overlay_img = 0;
        stride = 1;
    end
    img = EDTI.LoadNifti(nii_file);
    img = img.img;
    img = single(img);
    img = img/prctile(img(:),99);
    if(ax_or_coronal == 1)
      % Capture the plot as an image 
      if(ndims(img) == 3)
          for z=1:stride:size(img,3)
              imagesc(img(:,:,z),[0 1]);
              axis image;
              axis off;
              colormap(cmap)
              frame = getframe(gcf); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256); 
              % Write to the GIF File 
              if z == 1 
                  imwrite(imind,cm,img_name,'gif', 'Loopcount',inf); 
              else 
                  imwrite(imind,cm,img_name,'gif','WriteMode','append','DelayTime',0.1); 
              end  
          end
      elseif(ndims(img) == 4)
          for z=1:stride:size(img,4)
              imagesc(img(:,:,round(end/2),z),[0 1]);
              axis image;
              axis off;
              colormap(cmap)
              frame = getframe(gcf); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256); 
              % Write to the GIF File 
              if z == 1 
                  imwrite(imind,cm,img_name,'gif', 'Loopcount',inf); 
              else 
                  imwrite(imind,cm,img_name,'gif','WriteMode','append','DelayTime',0.1); 
              end           
          end
      end
    end
    
    [out_folder,out_name] = fileparts(img_name);
    new_nii_name = fullfile(nifti_folder,[out_name '.nii']);
    copyfile(nii_file,new_nii_name);
    if(overlay_img == 0)
        html_code = ['<a href="papaya.html?image_uri=' urlencode(strrep(new_nii_name,'\','/')) '" '... 
                img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    else
        html_code = ['<a href="papaya.html?overlay_uri=' urlencode(strrep(new_nii_name,'\','/')) ...
            '&image_uri=' urlencode(strrep(overlay_img,'\','/')) '" '...
                img_name '"><img src="' img_name ...
                '" href="' img_name '" '...
                ' style="width:' num2str(width_in_px) 'px" /></a>'];
    end
end