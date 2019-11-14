%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function nii_splitter(first_file_path,shells,out_prefix)
% This function to split by b-val 
% Alberto De Luca


decomp = strsplit(first_file_path,'.');
first_file = decomp{1};
img_one = [first_file '.nii.gz']; 
bvecs1 = load([first_file '.bvec']);
bvals1 = load([first_file '.bval']);

data = load_untouch_nii(first_file_path);
[sx,sy,sz,ss] = size(data.img);

B0s = find(bvals1<11);
for j=1:length(shells)
    disp(['---Splitting shell ' num2str(shells(j))]);
    shell_index = find(abs(shells(j)-bvals1)<110);
    out_shell = zeros(sx,sy,sz,length(shell_index)+length(B0s));
    out_shell(:,:,:,1:length(B0s)) = data.img(:,:,:,B0s);
    out_shell(:,:,:,length(B0s)+1:end) = data.img(:,:,:,shell_index);
    save_file.hdr = data.hdr;
    save_file.hdr.dime.dim(5) = size(out_shell,4);
    save_file.img = out_shell;
    try
        save_nii(save_file,[first_file '_' num2str(shells(j)) '.nii.gz']);
    catch err
        disp('Orig file header misses some fields, creating a new one.');
        save_nii(make_nii(out_shell,save_file.hdr.dime.pixdim(2:4)),[out_prefix '_' num2str(shells(j)) '.nii.gz']);
    end
    bvals_new = bvals1([B0s shell_index]);
    bvecs_new = bvecs1([B0s shell_index],:);
    save([out_prefix '_' num2str(shells(j)) '.bval'],'bvals_new','-ascii');
    save([out_prefix '_' num2str(shells(j)) '.bvec'],'bvecs_new','-ascii');
    disp('---Done---');
end
disp('All done');
