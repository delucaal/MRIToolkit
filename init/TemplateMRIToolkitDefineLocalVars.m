%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function MRIToolkitDefineLocalVars()
    global MRIToolkit;
    MRIToolkit.Elastix.Location = 'The location of Elastix for your OS'; % Adjust per OS
    elastix_cmd = dir(fullfile(MRIToolkit.Elastix.Location,'*elastix*'));
    MRIToolkit.Elastix.ElastixCMD = fullfile(MRIToolkit.Elastix.Location,elastix_cmd.name);
    transformix_cmd = dir(fullfile(MRIToolkit.Elastix.Location,'*transformix*'));
    MRIToolkit.Elastix.TransformixCMD = fullfile(MRIToolkit.Elastix.Location,transformix_cmd.name);
	MRIToolkit.dcm2niix = 'The location of the dcm2niix executble (it is also part of MRIcroGL)';
	MRIToolkit.JSONio = 'location of the JSONio library, available at https://github.com/gllmflndn/JSONio';
    MRIToolkit.spm_path = 'Path to SPM + CAT12';
end