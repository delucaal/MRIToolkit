%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%
function MRIToolkitDefineLocalVars()
    global MRIKaleidoscope;
    MRIKaleidoscope.Elastix.Location = '~/Desktop/M_Code_ExploreDTI_v4.8.6/Source/MD_cor_E/macOSX64'; % Adjust per OS
    elastix_cmd = dir(fullfile(MRIKaleidoscope.Elastix.Location,'*elastix*'));
    MRIKaleidoscope.Elastix.ElastixCMD = fullfile(MRIKaleidoscope.Elastix.Location,elastix_cmd.name);
    transformix_cmd = dir(fullfile(MRIKaleidoscope.Elastix.Location,'*transformix*'));
    MRIKaleidoscope.Elastix.TransformixCMD = fullfile(MRIKaleidoscope.Elastix.Location,transformix_cmd.name);
end