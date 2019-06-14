%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of CC BY-NC-ND (https://creativecommons.org/licenses) %%%
% A. De Luca - UMC Utrecht - alberto@isi.uu.nl
% This script performs a T1 fit with inversion recovery data
% This code uses the RelaxometryMethods and OptimizationMethods libraries
% We here solve the problem as: (1-2*exp(-params.TI/T1)+exp(-params.TR./T1))
% with a grid-search approach

clear

% Each file corresponds to a different inversion time
sequences_files = {
    'AD09_WIP_qT1_IR_50_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_100_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_200_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_400_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_800_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_1600_20190606165250_1_2_SSFP.nii',...
    'AD09_WIP_qT1_IR_3200_20190606165250_1_2_SSFP.nii'
};

% The first number is the inversion time
properties = {
        {'TI',0.050,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',0.100,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',0.200,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',0.400,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',0.800,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',1.600,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
        {'TI',3.200,'TR',5,'TE',0.0038,'AcqType',MRIacq.ACQTYPE_IR,'AcqMode',MRIacq.ACQMODE_3D},...
};

ir_acq = MRIacq.load_volumes_with_properties(sequences_files,properties);
ir_acq = ir_acq.sort_by('TI','ascend');
ir_acq = ir_acq.normalize_4d();
ir_acq = ir_acq.vectorize_data();

%%

% Minimum - Maximum allowed T1s during the grid search
T1_min_max = [0.1 3];
% Number of points to create a linear range between min-max. Note: the
% solution is discretized.
T1_points = 5000;

params = ir_acq.properties_struct();

tic

[T1_map,RSS] = RelaxometryMethods.IR_T1_fit_ND(T1_min_max,T1_points,ir_acq.img,params);

T1_map = T1_map(:,2);
msize = ir_acq.matrix_size();
T1_map = reshape(T1_map,msize(1),msize(2),msize(3));
RSS = reshape(RSS,msize(1),msize(2),msize(3));

toc

save('IR_T1_map','T1_map')