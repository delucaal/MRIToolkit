%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



% A. De Luca - UMC Utrecht - alberto@isi.uu.nl
% This script performs a T1 fit with inversion recovery data
% This code uses the RelaxometryMethods and OptimizationMethods libraries
% We here solve the problem as a linear monoexponential 
clear

% %%% Example if your multiple echoes are stored in multiple .nii files
% % Each file corresponds to a different echo time
% sequences_files = {
%     'AD09_WIP_qT2-MESE_20190606165250_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e2_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e3_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e4_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e5_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e6_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e7_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e8_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e9_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e10_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e11_2_SSFP.nii',...
%     'AD09_WIP_qT2-MESE_20190606165250_e12_2_SSFP.nii'
%     };
% 
% % The first number is the inversion time, the second the repetition time, the third the echo-time
% properties = {
%         {'TI',-1,'TR',2,'TE',0.02,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.04,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.06,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.08,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.10,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.12,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.14,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.16,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.18,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.20,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.22,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
%         {'TI',-1,'TR',2,'TE',0.24,'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS},...
% };
% 
% se_acq = MRIacq.load_volumes_with_properties(sequences_files,properties);
%%% End example

% Now let's deal with data in a single matrix

my_mat_file = load('t2_mapping.mat');

TR = zeros(size(my_mat_file.t2,4),1); % TR does not matter for this fit - just a placeholder
TI = zeros(size(my_mat_file.t2,4),1); % Same for TI

se_acq = MRIacq('img',single(my_mat_file.t2),'TE',my_mat_file.echo_time,...
    'AcqType',MRIacq.ACQTYPE_SE,'AcqMode',MRIacq.ACQMODE_2DMS,'TR',TR,'TI',TI);

se_acq = se_acq.sort_by('TE','ascend');
se_acq = se_acq.discard_volumes([1 2 3]);% discard the first three echoes due to stimulated echoes
se_acq = se_acq.vectorize_data();

params = se_acq.properties_struct();

tic

[r2_map,cost,fit] = RelaxometryMethods.SE_R2_fit_ND(se_acq.img,params);

msize = se_acq.matrix_size();
t2_map = 1./r2_map(:,2);
t2_map(~isfinite(t2_map) | t2_map > 10 | t2_map < 0) = 0;
t2_map = reshape(t2_map,msize(1),msize(2),msize(3));

cost = reshape(cost,msize(1),msize(2),msize(3));
fit = reshape(fit,msize(1),msize(2),msize(3),size(fit,2));

se_acq = se_acq.reshape_data_to_original();

toc

save('Monoexp_MESE_T2_map','t2_map');