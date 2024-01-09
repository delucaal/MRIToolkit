%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

function DW_Elastix_Register(moving,fixed,regparlist,output_dir,mMask,fMask,outfilename)
    global MRIToolkit;
    if(isfield(MRIToolkit,'Elastix'))
        elastix_path = MRIToolkit.Elastix.Location;
        elastix_cmd = MRIToolkit.Elastix.ElastixCMD;
        if(~ispc)
            if(~ismac)
                elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' elastix_cmd];
            else
                elastix_cmd = ['DYLD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'elastix_Mac_64'];
            end
        end
    else
        elastix_path = '/Users/alb/Desktop/M_Code_ExploreDTI_v4.8.6/Source/MD_cor_E/macOSX64/';
        if(~ismac)
            elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'elastix_Mac_64'];
        else
            elastix_cmd = ['DYLD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'elastix_Mac_64'];
        end
    end
    sentence = [elastix_cmd ' -m ' moving ' -f ' fixed ' -out ' output_dir];
    if(exist('mMask','var') > 0 && ~isempty(mMask))
       sentence = [sentence ' -mMask ' mMask];
    end
    if(exist('fMask','var') > 0 && ~isempty(fMask))
       sentence = [sentence ' -fMask ' fMask];
    end
    for j=1:length(regparlist)
       sentence = [sentence ' -p ' regparlist{j}]; 
    end
    if(~exist(output_dir,'dir'))
        mkdir(output_dir)
    end
    system(sentence);
    if(exist('outfilename','var') > 0)
        result = dir([output_dir '/result*.nii.gz']);
        if(isempty(result))
           tp = dir([output_dir '/TransformParameters.*.txt']);
           DW_Elastix_Transform(moving,outfilename,[output_dir '/' tp(end).name]); 
        else
           if(length(result) == 1)  
               try
                  copyfile(fullfile(output_dir,result.name),outfilename);
               catch
                   system(['cp ' fullfile(output_dir,result.name) ' ' outfilename]);
               end
           else
               copyfile(fullfile(output_dir,result(end).name),outfilename);
               system(['cp ' fullfile(output_dir,result(end).name) ' ' outfilename]);
           end
        end
    end
end