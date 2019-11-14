%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function MRIToolkitDefineLocalVars()
    global MRIToolkit;
    MRIToolkit.Elastix.Location = '~/Desktop/M_Code_ExploreDTI_v4.8.6/Source/MD_cor_E/macOSX64'; % Adjust per OS
    elastix_cmd = dir(fullfile(MRIToolkit.Elastix.Location,'*elastix*'));
    MRIToolkit.Elastix.ElastixCMD = fullfile(MRIToolkit.Elastix.Location,elastix_cmd.name);
    transformix_cmd = dir(fullfile(MRIToolkit.Elastix.Location,'*transformix*'));
    MRIToolkit.Elastix.TransformixCMD = fullfile(MRIToolkit.Elastix.Location,transformix_cmd.name);
end