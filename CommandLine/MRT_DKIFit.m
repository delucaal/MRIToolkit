function MRT_DKIFit(varargin)
    disp(varargin);
    path = varargin{1};
    files = dir(path);
    for file_id=1:length(files)
       disp(['Processing file ' num2str(file_id) ' of ' num2str(length(files))]);
       disp(files(file_id).name);
       
       basename = fullfile(files(file_id).folder,files(file_id).name);
       if(contains(basename,'.nii.gz'))
           [fp,fn] = fileparts(basename(1:end-3));
       else
           [fp,fn] = fileparts(basename);
       end
       basename = fullfile(fp,fn);
       
       EDTI.b_Matrix_from_bval_bvec([basename '.bval'],[basename '.txt']);
       EDTI.ReplaceInvalidKurtosisPoints(true);
% %        % TEST
% %        EDTI.SelectVolumesWithBvals('nii_file',[basename '.nii'],'output',fullfile(tempdir,'OnlyB2500.nii'),...
% %            'bvals',[0 2500]);
% %        EDTI.b_Matrix_from_bval_bvec(fullfile(tempdir,'OnlyB2500.bval'),fullfile(tempdir,'OnlyB2500.txt'));
% %        EDTI.PerformDTI_DKIFit('nii_file',fullfile(tempdir,'OnlyB2500.nii'),'grad_perm',2,'grad_flip',1,'dki',0,...
% %            'fit_mode','ols');
% %        copyfile(fullfile(tempdir,'OnlyB2500.mat'),[basename '.mat']);
% %        % TEST
       EDTI.PerformDTI_DKIFit('nii_file',[basename '.nii'],'grad_perm',2,'grad_flip',1,'dki',1,...
           'fit_mode','wls'...
           ,'dki_constraints','[0 Inf]');
       EDTI.MatMetrics2Nii([basename '.mat'],1);
       
    end    
end
