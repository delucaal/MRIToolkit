%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

function DW_Elastix_Transform(filein,fileout,transform_parameters)
    transformix_cmd = ImageRegistrations.BasicTransformixCall();
    mkdir('temp_transformix');  
    if(contains(filein,'.nii') < 1)
        header = load_untouch_header_only([filein '.nii']);
    else
        header = load_untouch_header_only(filein);
    end
    if(header.dime.dim(1) == 3)
        disp('3D file');
        sentence = [transformix_cmd ' -in ' filein ' -out temp_transformix/ -tp ' transform_parameters];
        system(sentence);
%         gzip('temp_transformix/*.nii');
        result_file = dir('temp_transformix/result*.nii*');
        if(strcmp(result_file.name(end),'z'))
           if(~strcmp(fileout(end-1:end),'gz'))
               fileout = [fileout '.gz'];
           end
        else
            if(strcmp(fileout(end-1:end),'gz'))
                fileout = fileout(1:end-3);
            end
        end
        copyfile(['temp_transformix/' result_file.name],fileout);
%         system(['cp temp_transformix/result.nii.gz ' fileout]);
    elseif(header.dime.dim(1) == 4)
        disp('4D file');
        data = load_untouch_nii(filein);
        temp_data = data;
        temp_data.hdr.dime.dim(1) = 3;
        temp_data.hdr.dime.dim(5) = 1;
        system('mkdir temp_transformix_4d');
        for vol_id=1:header.dime.dim(5)
            temp_data.img = data.img(:,:,:,vol_id);
            save_untouch_nii(temp_data,'./temp_transformix_4d/tmp_elastix_transform.nii.gz');
            out_name = fullfile('temp_transformix_4d',sprintf('reg_tmp_reg_elastix_%.3d.nii',vol_id));
            DW_Elastix_Transform('./temp_transformix_4d/tmp_elastix_transform.nii.gz',out_name,transform_parameters);
            J = load_untouch_nii(out_name);
            if(vol_id == 1)
                nfile = zeros(J.hdr.dime.dim(2:5));
            end
            nfile(:,:,:,vol_id) = J.img;
        end
        data.img = nfile;
        data.hdr.dime.dim(2:5) = size(data.img);
        save_untouch_nii(data,fileout);
%         system(['fslmerge -t ' fileout ' ./temp_transformix_4d/reg_tmp_reg_elastix*.gz']);
%         system('rm -rf temp_transformix_4d');
    end
    try
        delete('temp_transformix/*');
        rmdir('temp_transformix');
    catch
    end
end