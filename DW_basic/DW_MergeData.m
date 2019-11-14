%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function [new_data,parts] = DW_MergeData(data1,data2)
    if(ischar(data1))
        if(exist([data1 '.nii'],'file') > 0)
            data1 = DW_LoadData([data1 '.nii'],[data1 '.bvec'],[data1 '.bval']);
        elseif(exist([data1 '.nii.gz'],'file') > 0)
            data1 = DW_LoadData([data1 '.nii.gz'],[data1 '.bvec'],[data1 '.bval']);
        else
            error('Cannot find file 1');
        end
    end

    if(ischar(data2))
        if(exist([data2 '.nii'],'file') > 0)
            data2 = DW_LoadData([data2 '.nii'],[data2 '.bvec'],[data2 '.bval']);
        elseif(exist([data1 '.nii.gz'],'file') > 0)
            data2 = DW_LoadData([data2 '.nii.gz'],[data2 '.bvec'],[data2 '.bval']);
        else
            error('Cannot find file 2');
        end
    end
    
    parts = zeros(data1.hdr.dime.dim(5)+data2.hdr.dime.dim(5),1);
    parts(1:data1.hdr.dime.dim(5)) = 1;
    parts(1+data1.hdr.dime.dim(5):end) = 2;
    
    new_data = data1;
    new_data.hdr.dime.dim(5) = data1.hdr.dime.dim(5)+data2.hdr.dime.dim(5);
    new_data.img = zeros(new_data.hdr.dime.dim(2:5));
    new_data.img(:,:,1:data1.hdr.dime.dim(4),1:data1.hdr.dime.dim(5)) = data1.img;
    new_data.img(:,:,1:data2.hdr.dime.dim(4),1+data1.hdr.dime.dim(5):end) = data2.img;
    new_data.bvals = [data1.bvals;data2.bvals];
    new_data.bvecs = [data1.bvecs;data2.bvecs];
    if(isfield(data1,'mask'))
        new_data.mask = data1.mask;
    end
end