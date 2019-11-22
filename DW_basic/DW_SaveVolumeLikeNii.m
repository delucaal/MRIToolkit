%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



function DW_SaveVolumeLikeNii(input_vol,ref_nii_struct,out_prefix,niigz)
    if(nargin < 4)
        niigz = 1;
    end
    out_vol.img = input_vol;
    out_vol.hdr = ref_nii_struct.hdr;
    if(isfield(ref_nii_struct,'untouch'))
        out_vol.untouch = 1;
    end
    out_vol.hdr.dime.dim(1) = ndims(input_vol);
    if(ndims(input_vol) < 4)
        out_vol.hdr.dime.dim(5) = 1;
    end
    out_vol.hdr.dime.dim(2:1+out_vol.hdr.dime.dim(1)) = size(input_vol);
    
    if(isinteger(input_vol))
        out_vol.hdr.dime.bitpix = 16;
        out_vol.hdr.dime.datatype = 512;
    else
        out_vol.hdr.dime.bitpix = 32;
        out_vol.hdr.dime.datatype = 64;
    end
    out_vol.hdr.dime.scl_slope = 1;
    out_vol.hdr.dime.scl_inter = 0;
    if(niigz == 1)
        save_untouch_nii(out_vol,[out_prefix '.nii.gz']);
    else
        save_untouch_nii(out_vol,[out_prefix '.nii']);
    end
end