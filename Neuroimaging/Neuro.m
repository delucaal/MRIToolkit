%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% This class executes common automatic pipelines and structures their
% output in an MRIToolkit friendly way

classdef Neuro < handle
   methods (Static)
       
       function CAT12Pipeline(varargin)
            global MRIToolkit
            if(isempty(MRIToolkit.spm_path))
               error('No SPM path was specified. Please add it to the MRIToolkitDefineLocalVars.m file');               
            end
            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            output = GiveValueForName(coptions,'output');
            if(~contains(nii_file,'.gz') && exist(nii_file,'file') < 1)
                nii_file = [nii_file '.gz'];
            end
            [fp,fn,ext] = fileparts(nii_file);
            if(isempty(output))
                output = fullfile(fp,[fn '_CAT12']);
            end
            
%             if(contains(output,pwd) == false)
                mkdir(output);
                copyfile(nii_file,fullfile(output,[fn ext]));
%             end
            force_wmh = GiveValueForName(coptions,'force_wmh');
            if(isempty(force_wmh))
                force_wmh = 0;
            end
            
            flair_file = GiveValueForName(coptions,'flair_file');
            if(isempty(flair_file))
                flair_file = [];
            end            
            
            ExecuteCat12_T1_FLAIR(fullfile(output,[fn ext]),flair_file,0,0,force_wmh);
       end
       
       function CAT12ApplyDeformation(varargin)
            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            field_file = GiveValueForName(coptions,'field_file');
            if(isempty(field_file))
                error('Missing mandatory argument field_file');
            end
            
            matlabbatch{1}.spm.tools.cat.tools.defs2.field = {field_file};
            matlabbatch{1}.spm.tools.cat.tools.defs2.images = {{nii_file}};
            matlabbatch{1}.spm.tools.cat.tools.defs2.interp = 0;
            matlabbatch{1}.spm.tools.cat.tools.defs2.modulate = 0; 
            
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            
            [fp,fn,ext] = fileparts(nii_file);
            EDTI.FlipPermuteSpatialDimensions('nii_file',fullfile(fp,['w' fn '.nii']),'output',fullfile(fp,['w' fn '_FP.nii']),'flip',[0 1 0])
       end
       
       function SPMBiasFieldCorrection(varargin)
            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            
            ExecuteBiasFieldCorrection(nii_file);

       end
       
       function MakeAnimatedGIFOfVolume(varargin)
            coptions = varargin;
            gif_file = GiveValueForName(coptions,'gif_file');
            if(isempty(gif_file))
                error('Missing mandatory argument gif_file');
            end
            interactive = GiveValueForName(coptions,'interactive');
            if(isempty(interactive))
                interactive = 1;
            end
            vol_in = GiveValueForName(coptions,'vol');
            if(isempty(vol_in))
                error('Missing mandatory argument vol');
            end 
            overlay_in = GiveValueForName(coptions,'overlay');
            if(isempty(overlay_in))
                overlay_in = [];
            end
            ref_ax = GiveValueForName(coptions,'axes');
            if(isempty(overlay_in))
                ref_ax = 1;
            end           
            
            if(ref_ax == 1)
            elseif(ref_ax == 2)
                vol_in = permute(vol_in,[3 2 1]);
                vol_in = (flip(vol_in,1));
                overlay_in = permute(overlay_in,[3 2 1]);
                overlay_in = (flip(overlay_in,1));               
            elseif(ref_ax == 3)
                vol_in = permute(vol_in,[3 1 2]);
                vol_in = (flip(vol_in,1));
                overlay_in = permute(overlay_in,[3 1 2]);
                overlay_in = (flip(overlay_in,1));               
            end
            
            f=figure;
            filename = gif_file;
            set(gcf,'color','k','inverthardcopy','off')
            if(interactive == 0)
               set(gcf,'visible','off'); 
            end
            ax1=axes;
            ax2=axes;
            for iz=1:size(vol_in,3)
                imagesc(ax1,vol_in(:,:,iz));
                ax1.CLim = [0 max(vol_in(:))];
                colormap(ax1,'gray');
                axis(ax1,'image','off');
                if(~isempty(overlay_in))
                    h=imagesc(ax2,overlay_in(:,:,iz));
                    colormap(ax2,'hot');
                    ax2.CLim = [0 max(overlay_in(:))];
                    set(h,'AlphaData',overlay_in(:,:,iz) ~= 0)
                end
                axis(ax2,'image','off');
                if(interactive == 1)
                    pause(0.1)
                end
                frame = getframe(f); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                % Write to the GIF File 
                if iz == 1 
                  imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
                else 
                  imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
                end 

            end

            close(f)
 
       end

       % This function is meant to perform brain extraction on T1/T2/FLAIR
       % images by performing a non-linear registration of the MNI template
       function BrainExtractionWithRegistration(varargin)
            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            
            while(true)
                temp_direc = fullfile(tempdir,['MRT_' num2str(randi(1e5))]);
                if(exist(temp_direc,'dir') < 1)
                    mkdir(temp_direc)
                    break
                end
            end
            
            DW_Elastix_Register(which('MNI152_T1_2mm.nii.gz'),nii_file,{which('parrig_EPI.txt'), which('paraff_EPI.txt'), which('parbsp_EPI.txt')},...
                temp_direc,[],[],fullfile(temp_direc,'final_trafo.nii'));

            DW_Elastix_Transform(which('MNI152_T1_2mm_brain.nii.gz'),fullfile(temp_direc,'final_trafo.nii'),...
                fullfile(temp_direc,'TransformParameters.2.txt'));

            mask = MRTQuant.LoadNifti(fullfile(temp_direc,'final_trafo.nii'));
            in_file = MRTQuant.LoadNifti(nii_file);
            in_file.img = single(in_file.img) .* single(mask.img > 0.05*max(mask.img(:)));
            
            MRTQuant.WriteNifti(in_file,output);
            rmdir(temp_direc,'s');
       end

       % Reconstruct structural connectivity given whole brain tractography
       % and corresponding brain parcellations in native space.
       function StructuralConnectivityAnalysis(varargin)
            coptions = varargin;
            mat_file = GiveValueForName(coptions,'mat_file');
            if(isempty(mat_file))
                error('Missing mandatory argument mat_file');
            end
            tract_file = GiveValueForName(coptions,'tract_file');
            if(isempty(tract_file))
                error('Missing mandatory argument tract_file');
            end
            label_file = GiveValueForName(coptions,'label_file');
            if(isempty(label_file))
                error('Missing mandatory argument label_file');
            end
            labels_txt = GiveValueForName(coptions,'labels_txt');
            if(isempty(labels_txt))
                error('Missing mandatory argument labels_txt');
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            if(exist(output,'dir') < 1)
                mkdir(output);
            end
            EDTI_Library.E_DTI_Network_analysis_from_ROI_L({mat_file},{tract_file},{label_file},labels_txt,1,...
                'pass',output,0);
        
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


function ExecuteCat12_T1_FLAIR(t1_file,flair_file,ncores,showreport,force_wmh)
global MRIToolkit;
% global cat;
spm_path = MRIToolkit.spm_path;

disp(['ExecuteCat12_T1_FLAIR T1:' t1_file ' FLAIR:' flair_file]);

if(isempty(which('spm')))
    addpath(genpath(spm_path));
end
cat12('expert')

if(nargin < 2 || showreport > 0)
    showreport = 2;
end

if(nargin < 3)
    ncores = 0;
    showreport = 0;
    force_wmh = 0;
end

matlabbatch{1}.spm.tools.cat.estwrite.data = {t1_file};
if(~isempty(flair_file) || force_wmh == 1)
    matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {flair_file};
end
matlabbatch{1}.spm.tools.cat.estwrite.nproc = ncores;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {fullfile(spm_path,'tpm','TPM.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 3;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.shootingtpm = {fullfile(spm_path,'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_0_GS.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb = 12;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP = 22;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1; % Changed on 31-10-2024 before processing HARMY
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
if(~isempty(flair_file) || force_wmh == 1)
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 1;
else
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
end
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

% save('test_job','matlabbatch')

spm('defaults', 'FMRI');
try
    spm_jobman('run', matlabbatch);
catch err
   disp('Error in CAT12');
   global cat
   disp(cat)
   disp(err.message);
end
end

function ExecuteBiasFieldCorrection(input_file)
    global MRIToolkit;
    spm_path = MRIToolkit.spm_path;
    
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {[input_file ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm_path,'/tpm/TPM.nii,1')};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm_path,'/tpm/TPM.nii,2')};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm_path,'/tpm/TPM.nii,3')};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm_path,'/tpm/TPM.nii,4')};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm_path,'/tpm/TPM.nii,5')};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm_path,'/tpm/TPM.nii,6')};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
end
