%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function data = DW_LoadData(data_name,bvec_name,bval_name,mask_name,load_nii_instead,noise_map,round_small_bs_to_zero)
if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
    data = load_nii(data_name);
    disp('Loading with load_nii');
else
    data = load_untouch_nii(data_name);
    disp('Loading with load_untouch_nii');
end

data.bvecs = load(bvec_name);
data.bvals = load(bval_name);
data.bvecs = data.bvecs';
data.bvals = data.bvals';
data.img(data.img<0) = 0;
if(data.hdr.dime.scl_slope == 0)
    data.hdr.dime.scl_slope = 1;
end
data.img = single(data.img)*data.hdr.dime.scl_slope+data.hdr.dime.scl_inter;

if(exist('round_small_bs_to_zero','var') > 0 && round_small_bs_to_zero > 0)
    disp('Rounding small bs to zero');
    data.bvals(data.bvals < 1) = 0;
end

if(exist('mask_name','var') > 0 && ~isempty(mask_name))
    mask = load_untouch_nii(mask_name);
    data.mask = mask.img;
else
    data.mask = data.img(:,:,:,find(data.bvals<1,1))>0;
end
if(size(data.img,4) > length(data.bvals))
%     WE HAVE NOISE MAP AT THE END
    disp('Noisemap found');
    data.noisemap = data.img(:,:,:,end);
    data.img = data.img(:,:,:,1:end-1);
    data.hdr.dime.dim(5) = size(data.img,4);
else
    if(exist('noise_map', 'var') > 0 && ~isempty(noise_map))
        disp('Loading noise map');
        try
        if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
            nmap = load_nii(noise_map);
        else
            nmap = load_untouch_nii(noise_map);
        end
        data.noisemap = nmap.img;
        catch
            disp('Error loading noisemap, please check ');
        end
    else
        %INDEED CHECK IF A NOISE MAP EXIST
        base_name = strsplit(data_name,'.');
        if(exist([base_name{1} '_noisemap.nii.gz'],'var') > 0)
           disp('Found a noisemap matching dataname, loading'); 
            if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
                nmap = load_nii([base_name{1} '_noisemap.nii.gz']);
            else
                nmap = load_untouch_nii([base_name{1} '_noisemap.nii.gz']);
            end
            data.noisemap = nmap.img;
        end
    end
end
end