% 1 Download data from https://surfdrive.surf.nl/files/index.php/s/kAYfZ0xugHKOC4t
% 2 Execute Example1 and generate some results

if(isempty(which('ReportMaker')))
    addpath('MRIToolkit');
    MRIToolkitInit;
end

RPM = ReportMaker();

RPM.AddDataFolder(pwd);
RPM.AddFolders('RawNii');
% RPM.AddFolder('derivatives/sub-Z_UT_424_ses-1');

% Add the T1
label = 'T1';
rel_path_to_file = 'anat';
search_string = '*_brain_FP_ds.nii*';
axial_coronal_sagittal = 1;
link_to_viewer = 0;
make_gif = 1;
multi_slice = 0;
cmap = colormap('gray');
overlay_on = 0;
stride = 5;
RPM.AddImage2Scan(label,rel_path_to_file,search_string,axial_coronal_sagittal,link_to_viewer,make_gif,multi_slice,cmap,overlay_on,stride);

% FA
label = 'FA';
rel_path_to_file = 'dwi';
search_string = '*trafo_FA.nii*';
axial_coronal_sagittal = 1;
link_to_viewer = 0;
make_gif = 1;
multi_slice = 0;
cmap = colormap('gray');
overlay_on = 1;
stride = 5;
RPM.AddImage2Scan(label,rel_path_to_file,search_string,axial_coronal_sagittal,link_to_viewer,make_gif,multi_slice,cmap,overlay_on,stride);

RPM.Process('TestReport')
