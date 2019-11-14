%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function DW_SaveData(data,out_name,niigz)
    if(nargin < 3)
        niigz = 1;
    end
    if(niigz == 1)
        ext = '.nii.gz';
    else
        ext = '.nii';
    end
    data.img = data.img/data.hdr.dime.scl_slope-data.hdr.dime.scl_inter;
    save_untouch_nii(data,[out_name ext]);
    if(isfield(data,'bvals'))
        bvals = data.bvals';
        bvecs = data.bvecs';
        save([out_name '.bvec'],'bvecs','-ascii');
        save([out_name '.bval'],'bvals','-ascii');
    end
    if(isfield(data,'mask') > 0)
       mask.hdr = data.hdr;
       mask.untouch = 1;
       mask.hdr.dime.dim(1) = 3;
       mask.hdr.dime.dim(5) = 1;
       mask.img = zeros(mask.hdr.dime.dim(2:4));
       mask.img = data.mask;
       save_untouch_nii(mask,[out_name '_mask' ext]);
    end
    if(isfield(data,'noisemap') > 0)
       mask.hdr = data.hdr;
       mask.untouch = 1;
       mask.hdr.dime.dim(1) = 3;
       mask.hdr.dime.dim(5) = 1;
       mask.img = zeros(mask.hdr.dime.dim(2:4));
       mask.img = data.noisemap;
       save_untouch_nii(mask,[out_name '_noisemap' ext]);
    end
end