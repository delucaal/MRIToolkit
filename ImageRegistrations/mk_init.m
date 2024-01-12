%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%% A. De Luca% Sub-toolbox specific initialization% 29/03/2018: creation - v1.0% 12/01/2024: v1.2 Added the ImageRegistrations classglobal MRIToolkit;MRIToolkit.image_registrations_version = 1.2;addpath(get_executed_file_path())addpath(fullfile(get_executed_file_path(),'elastix_parameters'))MRIToolkit.Elastix.Works = false;try    if(~isfield(MRIToolkit.Elastix,'Location') || exist(MRIToolkit.Elastix.Location,'dir') < 1)        if(exist(fullfile(userpath,'Elastix'),'dir') > 0)            MRIToolkit.Elastix.Location = fullfile(userpath,'Elastix');            MRIToolkit.Elastix.ElastixCMD = fullfile(userpath,'Elastix','elastix');            MRIToolkit.Elastix.TransformixCMD = fullfile(userpath,'Elastix','transformix');        else            MRIToolkit.Elastix.Location = '';        end    end    if(exist(MRIToolkit.Elastix.Location,'dir') > 0)        if(exist(MRIToolkit.Elastix.ElastixCMD,'file') > 0 && exist(MRIToolkit.Elastix.TransformixCMD,'file') > 0)            [a,b] = system(ImageRegistrations.BasicElastixCall());            [c,d] = system(ImageRegistrations.BasicTransformixCall());            if(contains(b,'elastix-usage') > 0 && contains(d,'transformix-usage') > 0)                MRIToolkit.Elastix.Works = true;            end        end    end    if(MRIToolkit.Elastix.Works == false)        warning('There seems to be a problem with the installation of Elastix. I can attempt fixing it by downloading elastix in your home folder. Call ImageRegistrations.DownloadElastix() to do so.');        warning('Motion correction and image registrations will not work until Elastix is properly configured.');    endcatchend