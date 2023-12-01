%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



% This class transforms the processing methods originally implemented in ExploreDTI
% into a library, made consistent into a class enforcing its consistency,
% ease of use and availability for scripting / command line tools without
% need for a GUI.
% Author: Alberto De Luca - alberto@isi.uu.nl - alberto.deluca.06@gmail.com
% Many of the methods here included were originally implemented by Alexander Leemans.
% First version 29/10/2019 as part of the class EDTI
% 24/10/2020: All ExploreDTI are adapted and moved to this class as static
% methods
classdef EDTI_Library < handle
    
    methods(Static)
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function DTI_files = E_DTI_SMECEPI_Get_input_DTI(par)
            
            if exist(par.DTI_f_in,'file')==2
                DTI_files{1} = par.DTI_f_in;
            elseif exist(par.DTI_f_in,'dir')==7
                files = EDTI_Library.E_DTI_Get_files_from_folder(par.DTI_f_in,'.mat');
                DTI_files = EDTI_Library.E_DTI_Get_DTI_files(files);
            else
                DTI_files=[];
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function files = E_DTI_Get_DTI_files(f)
            
            files = {};
            
            warning off all
            
            cn = 0;
            
            for i=1:length(f)
                
                try
                    load(f{i},'g','VDims','MDims','b','bval','NrB0')
                    if exist('g','var') && exist('NrB0','var') && exist('VDims','var') && exist('MDims','var') && exist('b','var') && exist('bval','var')
                        cn=cn+1;
                        files{cn} = f{i};
                        clear g VDims MDims b bval NrB0
                    end
                catch
                end
                
            end
            
            warning on all
        end
        
        % From ExploreDTI: helper function for masking
        function M = E_DTI_mean_DWI(DWI, NrB0)
            
            if ~iscell(DWI)
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            end
            
            pp = size(DWI{1});
            M = repmat(double(0),pp);
            
            for i=NrB0+1:length(DWI)
                M = M + double(DWI{i});
            end
            
            M = M/(length(DWI)-NrB0);
            
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function suc = E_DTI_do_initial_check_reg_tra_files(E_path,T_path)
            
            % disp(['E_path: ' E_path])
            
            suc = 1;
            
            if ispc
                sys_c = ['"' E_path '"'];
                [r,st]=system(sys_c);
            else
                try
                    [dir_e,el] = fileparts(E_path);
                    %     sys_c = ['./' el];
                    %     tedi = cd(dir_e);
                    %     [r,st]=system(sys_c);
                    %     cd(tedi)
                    if(ismac)
                        sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' fullfile(dir_e,el) '"'];
                    else
                        sys_c = ['LD_LIBRARY_PATH=' dir_e ' ' fullfile(dir_e,el)];
                    end
                    [r,st] = system(sys_c);
                catch me
                    disp('Debugging info for Alexander (please, email to: Alexander@isi.uu.nl)')
                    disp(['E_path: ' E_path])
                    disp(['el: ' el])
                    disp(['dir_e: ' dir_e])
                    disp(me.message)
                    suc=0;
                    return;
                end
            end
            if r~=0
                suc = 0;
                disp('Error related to file:')
                disp(E_path)
                disp('Error message:')
                disp(st);
            end
            
            if ispc
                sys_c = ['"' T_path '"'];
                [r,st]=system(sys_c);
            else
                [dir_e,tra] = fileparts(T_path);
                %     sys_c = ['./' tra];
                %     tedi = cd(dir_e);
                %     [r,st]=system(sys_c);
                %     cd(tedi)
                if(ismac)
                    sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' fullfile(dir_e,tra) '"'];
                else
                    sys_c = ['LD_LIBRARY_PATH="' dir_e '" "' fullfile(dir_e,tra) '"'];
                end
                [r,st] = system(sys_c);
            end
            if r~=0
                suc = 0;
                disp('Error related to file:')
                disp(T_path)
                disp('Error message:')
                disp(st);
            end
            
            if r~=0
                disp('For a solution, check the forum...')
                disp(' ')
            end
            
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function E_DTI_SMECEPI_Single(f_in,par)
            
            [FOL_FI,F,ext] = fileparts(f_in);
            
            tic;
            [suc, for_trafo] = EDTI_Library.E_DTI_SMECEPI_Single_only_EC(f_in,par);
            
            if suc==0
                return;
            end
            
            if par.R2D.type~=0
                
                EDTI_Library.E_DTI_SMECEPI_Single_only_EPI(f_in,for_trafo);
                
            end
            
            t=toc;
            if t<3600
                m = t/60;
                disp(['Computation time for ''' [F ext] ''' was ' num2str(m) ' minutes.'])
            elseif t>3600 && t<(3600*24)
                h = t/3600;
                disp(['Computation time for ''' [F ext] ''' was ' num2str(h) ' hours.'])
            else
                d = t/(3600*24);
                disp(['Computation time for ''' [F ext] ''' was ' num2str(d) ' days.'])
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % This actually performs the motion and eddy currents correction
        function [suc, for_trafo] = E_DTI_SMECEPI_Single_only_EC(f_in,par)
            
            for_trafo = struct;
            
            [FOL_FI,F] = fileparts(f_in);
            
            dir_temp = [par.temp_folder filesep 'Temp_' F];
            
            for_trafo.dir_temp = dir_temp;
            
            warning off all
            suc = mkdir(dir_temp);
            warning on all
            
            if suc==0
                disp('Error: could not write files/folders!')
                disp('Change permission properties of temp folder:')
                disp(par.temp_folder)
                return;
            end
            
            if ispc
                [suc, mess] = fileattrib(dir_temp,'+w -a -s','','s');
            else
                [suc, mess] = fileattrib(dir_temp,'+w +x','a','s');
            end
            
            if suc==0
                disp(['Problem with folder: ' dir_temp])
                disp(['Error message: ' mess])
                EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                return;
            end
            
            if par.DOF~=1
                par.Par_SMEC = [dir_temp filesep 'par_file_SMEC.txt'];
                
                suc = EDTI_Library.E_DTI_Par_MDC_aff_DTI(par);
                if suc==0
                    EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                    return;
                end
            end
            
            load(f_in,'VDims','DWI','FA')
            
            if par.DOF~=1
                fixed_SMEC_mask = single(~isnan(FA));
                
                if par.use_f_mask == 1
                    a=3;
                    se = zeros(a,a,a);
                    se((a+1)/2,(a+1)/2,(a+1)/2)=1;
                    se = smooth3(se,'gaussian',[a a a],1);
                    se = se>=se((a+1)/2,(a+1)/2,end);
                    fixed_SMEC_mask = single(imdilate(fixed_SMEC_mask,se)>0);
                else
                    fixed_SMEC_mask(:)=1;
                end
            end
            
            
            DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            
            for i=1:length(DWI)
                DWI{i} = single(DWI{i});
            end
            
            if par.DOF~=1
                par.fixed_SMEC_fn = [dir_temp filesep 'Fixed_SMEC.nii'];
                par.fixed_SMEC_fn_mask = [dir_temp filesep 'Fixed_SMEC_mask.nii'];
                
                fixed_SMEC = single(DWI{1});
                
                if par.Regul==1
                    fixed_SMEC = EDTI_Library.E_DTI_Anisotr_Gauss_smoothing(fixed_SMEC,VDims,2*max(VDims),[3 3 3]);
                end
                
                EDTI_Library.E_DTI_write_nifti_file(fixed_SMEC,VDims,par.fixed_SMEC_fn)
                EDTI_Library.E_DTI_write_nifti_file(fixed_SMEC_mask,VDims,par.fixed_SMEC_fn_mask)
                
                if exist(par.fixed_SMEC_fn,'file')~=2
                    suc = 0;
                    EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                    return;
                end
            end
            
            Le = length(DWI);
            
            f_moving = cell(Le,1);
            result = cell(Le,1);
            if par.DOF~=1
                DM_info = cell(Le,1);
            end
            VoxS = cell(Le,1);
            dir_temp_i = cell(Le,1);
            trafo_names = cell(Le,1);
            
            for i=1:Le
                
                dir_temp_i{i} = [dir_temp filesep 'Temp_' num2str(i)];
                warning off all
                suc = mkdir(dir_temp_i{i});
                warning on all
                
                if suc==0
                    disp('Error: could not write folder:')
                    disp(dir_temp_i{i})
                    EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                    return;
                end
                
                VoxS{i} = VDims;
                trafo_names{i} = [dir_temp_i{i} filesep 'TransformParameters.0.txt'];
                result{i} = [dir_temp_i{i} filesep 'result.nii'];
                f_moving{i} = [dir_temp_i{i} filesep 'moving_' num2str(i) '.nii'];
            end
            
            
            
            parfor i=1:Le
                EDTI_Library.E_DTI_output_TR_files(DWI{i},VoxS{i},f_moving{i},par)
            end
            
            if par.DOF~=1
                
                suc = EDTI_Library.E_DTI_SMECEPI_make_unity_transf(trafo_names{1},VoxS{1},size(DWI{1}));
                if suc==0
                    return;
                end
                
                parfor i=2:Le
                    EDTI_Library.E_DTI_do_the_SMEC_reg_step(f_moving{i},dir_temp_i{i},par)
                end
                
                if par.Regul==1
                    parfor i=2:Le
                        EDTI_Library.E_DTI_write_nifti_file(DWI{i},VoxS{i},f_moving{i});
                    end
                end
                
                parfor i=2:Le
                    EDTI_Library.E_DTI_do_the_SMEC_trafo_step(f_moving{i},dir_temp_i{i},trafo_names{i},par)
                end
                
                parfor i=2:Le
                    DWI{i} = EDTI_Library.E_DTI_read_nifti_file(result{i});
                end
                
                for i=1:length(DWI)
                    if isempty(DWI{i})
                        disp('Errors encountered...')
                        EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                        suc = 0;
                        return;
                    end
                end
                
                for i=1:length(DWI)
                    DWI{i} = single(DWI{i});
                    DWI{i}(DWI{i}==-1000)=nan;
                end
                
                for i=1:Le
                    
                    A = textread(trafo_names{i},'%s');
                    
                    DM_info{i}{1} = ...
                        [str2num(A{7})*(180/pi) str2num(A{6})*(180/pi) str2num(A{8})*(180/pi);
                        str2num(A{16}) str2num(A{15}) str2num(A{17}(1:end-1));
                        str2num(A{13}) str2num(A{12}) str2num(A{14});
                        str2num(A{10})  str2num(A{9}) str2num(A{11})];
                    
                    DM_info{i}{2} = EDTI_Library.E_DTI_Tra_Par_2_Tra_Mat(DM_info{i}{1});
                    
                    
                    DWI{i} = DWI{i}*det(DM_info{i}{2}(1:3,1:3));
                    
                end
                
                load(f_in,'b','NrB0','bval','info')
                
                b_old = b;
                
                for i=1:length(DWI)
                    World_Trafo = EDTI_Library.E_DTI_Tra_Par_2_Tra_Mat(DM_info{i}{1});
                    World_T = World_Trafo(1:3,1:3);
                    R = ((World_T*World_T')^(-1/2))*World_T;
                    B = [b_old(i,1) b_old(i,2)/2 b_old(i,3)/2;...
                        b_old(i,2)/2 b_old(i,4) b_old(i,5)/2;...
                        b_old(i,3)/2 b_old(i,5)/2 b_old(i,6)];
                    B_rot = R*B*R';
                    b(i,:) = [B_rot(1,1) 2*B_rot(1,2) 2*B_rot(1,3)...
                        B_rot(2,2) 2*B_rot(2,3) B_rot(3,3)];
                end
                
                if isnan(bval)
                    diff_model=2;
                else
                    diff_model=1;
                end
                
                [dummy, g] =  EDTI_Library.E_DTI_Get_Gradients_and_Bvalue(b, NrB0, diff_model);
                
                dummy = false(size(DWI{1}));
                for i=1:length(DWI)
                    dummy = or(dummy,isnan(DWI{i}));
                    DWI{i}(isnan(DWI{i}))=0;
                    DWI{i}(DWI{i}<0)=0;
                end
                
                if isempty(par.cust_mask.NS)
                    mask = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,par.mask_P.NS.NDWI,par.mask_P.NS.DWI,par.mask_P.NS.mfs);
                else
                    fn_cm = [FOL_FI filesep F par.cust_mask.NS];
                    [mask, VDimsm, suc] = EDTI_Library.E_DTI_read_nifti_file(fn_cm);
                    if suc==0
                        return;
                    end
                    mask = mask>0;
                end
                mask(dummy)=0;
                
                par_temp = par;
                par_temp.TE = par.TE.NS;
                
                if diff_model==1
                    [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par_temp,NrB0);
                elseif diff_model==2
                    [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par_temp,NrB0,VDims);
                end
                
                [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
                
                g = g(NrB0+1:end,:);
                
                if ~isreal(FA)
                    FA = real(FA);
                end
                if ~isreal(FEFA)
                    FEFA = real(FEFA);
                end
                if ~isreal(FE)
                    FE = real(FE);
                end
                if ~isreal(SE)
                    SE = real(SE);
                end
                if ~isreal(eigval)
                    eigval = real(eigval);
                end
                
                
                MDims = size(mask);
                
                f_out = [par.out_folder filesep F par.suff.NS];
                
                for_trafo.f_out = f_out;
                
                max_DWI = 0;
                for i=1:length(DWI)
                    max_DWI = max(max_DWI,max(DWI{i}(:)));
                end
                
                if max_DWI<=intmax('int16')
                    for i=1:length(DWI)
                        DWI{i} = round(DWI{i});
                        DWI{i} = int16(DWI{i});
                    end
                elseif max_DWI<=intmax('uint16')
                    for i=1:length(DWI)
                        DWI{i} = round(DWI{i});
                        DWI{i}(DWI{i}<0)=0;
                        DWI{i} = uint16(DWI{i});
                    end
                end
                
                try
                    save(f_out,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
                        'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','DM_info','par')
                catch me
                    suc = 0;
                    disp(me.message)
                    return;
                end
                
                if diff_model==2
                    save(f_out,'KT','-append')
                end
                
                if par.R2D.type==0
                    EDTI_Library.E_DTI_remove_temp_f(dir_temp);
                end
                
            else
                
                for_trafo.f_out=[];
                for i=1:Le
                    
                    suc = EDTI_Library.E_DTI_SMECEPI_make_unity_transf(trafo_names{i},VoxS{i},size(DWI{i}));
                    if suc==0
                        return;
                    end
                end
                
            end
            
            for_trafo.f_moving = f_moving;
            for_trafo.dir_temp_i = dir_temp_i;
            for_trafo.trafo_names = trafo_names;
            for_trafo.result_names = result;
            for_trafo.par = par;
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % This actually performs the EPI correction
        function E_DTI_SMECEPI_Single_only_EPI(f_in,for_trafo)
            
            [FOLD,FFO] = fileparts(f_in);
            
            par = for_trafo.par;
            
            [suc, FN_nii] = EDTI_Library.E_DTI_SMECEPI_check_data_stuff(f_in,par);
            
            if suc==0
                return;
            end
            
            fn_fixed = [for_trafo.dir_temp filesep 'Trafo_fixed.nii'];
            fn_fixed_mask = [for_trafo.dir_temp filesep 'Trafo_fixed_mask.nii'];
            [Fixed,VDims] = EDTI_Library.E_DTI_read_nifti_file(FN_nii);
            VDims_E = VDims;
            Fixed = single(Fixed);
            Fixed = 10000*(Fixed/max(Fixed(:)));
            EDTI_Library.E_DTI_write_nifti_file(Fixed,VDims,fn_fixed);
            
            mask = single(Fixed>0);
            if par.EPI.use_f_mask~=1
                mask(:)=1;
            end
            EDTI_Library.E_DTI_write_nifti_file(mask,VDims,fn_fixed_mask);
            
            fn_moving = [for_trafo.dir_temp filesep 'Trafo_moving.nii'];
            
            if ~isempty(for_trafo.f_out)
                LoadF = for_trafo.f_out;
            else
                LoadF = f_in;
            end
            
            if par.R2D.contrast==1
                load(LoadF,'DWI','VDims')
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
                Moving = single(DWI{1});
                clear DWI;
            elseif par.R2D.contrast==2
                load(LoadF,'FA','VDims')
                FA(isnan(FA))=0;
                FA = FA/sqrt(3);
                Moving = single(FA);
                clear FA;
            elseif par.R2D.contrast==3
                load(LoadF,'DWI','VDims','NrB0')
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
                Moving = single(E_DTI_mean_DWI(DWI, NrB0));
                clear DWI;
            end
            
            yp = prctile(Moving(:),99);
            Moving(Moving>yp)=yp;
            
            Moving = 10000*(Moving/max(Moving(:)));
            
            EDTI_Library.E_DTI_write_nifti_file(Moving,VDims,fn_moving);
            
            
            fn_par_trafo_rigid = [for_trafo.dir_temp filesep 'Par_file_Trafo_Rigid.txt'];
            
            par_EPI.Hist_bin = par.EPI.Hist_bin;
            par_EPI.Num_Resol = par.EPI.Num_Resol;
            par_EPI.Num_iter = par.EPI.Num_iter;
            par_EPI.Num_samp = par.EPI.Num_samp;
            par_EPI.Interpol = par.Interpol;
            par_EPI.Grid_Spacing = par.EPI.Grid_Spacing([2 1 3]);
            par_EPI.Par_FN = fn_par_trafo_rigid;
            par_EPI.Deriv_Scales = par.EPI.Deriv_Scales([2 1 3]);
            par_EPI.failsafe = par.R2D.failsafe;

            if(isfield(for_trafo.par,'use_NC') && for_trafo.par.use_NC == 1)
                suc = EDTI_Library.E_DTI_Par_Trafo_DTI_rigid_NC(par_EPI);
            else
                suc = EDTI_Library.E_DTI_Par_Trafo_DTI_rigid(par_EPI);
            end
            
            if suc==0
                EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                return;
            end
            
            if par.R2D.type~=3
                
                if ispc
                    
                    sys_c = ['"' par.E_path '"' ' -f ' '"' fn_fixed '"'...
                        ' -fMask ' '"' fn_fixed_mask '"' ...
                        ' -m ' '"' fn_moving '"'...
                        ' -out ' '"' for_trafo.dir_temp '"' ...
                        ' -p ' '"' fn_par_trafo_rigid '"'];
                    [r,st]=system(sys_c);
                    
                else
                    
                    [dir_e,el] = fileparts(par.E_path);
                    %         sys_c = ['./' el ' -f ' fn_fixed ...
                    %             ' -fMask ' fn_fixed_mask ...
                    %             ' -m ' fn_moving ...
                    %             ' -out ' for_trafo.dir_temp ...
                    %             ' -p ' fn_par_trafo_rigid];
                    if(ismac)
                        sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' fn_fixed ...
                            ' -fMask ' fn_fixed_mask ...
                            ' -m ' fn_moving ...
                            ' -out ' for_trafo.dir_temp ...
                            ' -p ' fn_par_trafo_rigid];
                    else
                        sys_c = ['LD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' fn_fixed ...
                            ' -fMask ' fn_fixed_mask ...
                            ' -m ' fn_moving ...
                            ' -out ' for_trafo.dir_temp ...
                            ' -p ' fn_par_trafo_rigid];
                    end
                    %         tedi = cd(dir_e);
                    [r,st]=system(sys_c);
                    %         cd(tedi)
                end
                
            else
                
                fn_par_trafo_non_rigid = [for_trafo.dir_temp filesep 'Par_file_Trafo_Non_Rigid.txt'];
                
                par_EPI.Par_FN = fn_par_trafo_non_rigid;
                if(isfield(for_trafo.par,'use_NC') && for_trafo.par.use_NC == 1)
                    suc = EDTI_Library.E_DTI_Par_Trafo_DTI_non_rigid_NC(par_EPI);
                else
                    suc = EDTI_Library.E_DTI_Par_Trafo_DTI_non_rigid(par_EPI);
                end
                
                if suc==0
                    EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                    return;
                end
                
                if ispc
                    
                    sys_c = ['"' par.E_path '"' ' -f ' '"' fn_fixed '"'...
                        ' -fMask ' '"' fn_fixed_mask '"' ...
                        ' -m ' '"' fn_moving '"'...
                        ' -out ' '"' for_trafo.dir_temp '"' ...
                        ' -p ' '"' fn_par_trafo_rigid '"' ...
                        ' -p ' '"' fn_par_trafo_non_rigid '"'];
                    [r,st]=system(sys_c);
                    
                else
                    
                    [dir_e,el] = fileparts(par.E_path);
                    %         sys_c = ['./' el ' -f ' fn_fixed ...
                    %             ' -fMask ' fn_fixed_mask ...
                    %             ' -m ' fn_moving ...
                    %             ' -out ' for_trafo.dir_temp ...
                    %             ' -p ' fn_par_trafo_rigid ...
                    %             ' -p ' fn_par_trafo_non_rigid];
                    if(ismac)
                        sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' fn_fixed ...
                            ' -fMask ' fn_fixed_mask ...
                            ' -m ' fn_moving ...
                            ' -out ' for_trafo.dir_temp ...
                            ' -p ' fn_par_trafo_rigid ...
                            ' -p ' fn_par_trafo_non_rigid];
                    else
                        sys_c = ['LD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' fn_fixed ...
                            ' -fMask ' fn_fixed_mask ...
                            ' -m ' fn_moving ...
                            ' -out ' for_trafo.dir_temp ...
                            ' -p ' fn_par_trafo_rigid ...
                            ' -p ' fn_par_trafo_non_rigid];
                    end
                    %         tedi = cd(dir_e);
                    [r,st]=system(sys_c);
                    %         cd(tedi)
                end
                
                
            end
            
            
            if r~=0
                disp('Error(s) during motion/distortion correction:')
                disp(' ')
                disp(st);
                disp(' ')
                disp('See the forum for a solution...')
                EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                return;
            end
            
            Trafo_rig_result = [for_trafo.dir_temp filesep 'TransformParameters.0.txt'];
            % Rig = EDTI_Library.E_DTI_get_Trafo_par_EPI_rigid_rotation_components(Trafo_rig_result);
            
            Q = textread(Trafo_rig_result,'%s');
            Rig{1} = str2num(Q{7});
            Rig{2} = str2num(Q{6});
            Rig{3} = str2num(Q{8});
            
            
            At = textread(Trafo_rig_result,'%s','delimiter','\n','bufsize',2^18);
            
            if par.R2D.type==3
                
                Trafo_nonrig_result = [for_trafo.dir_temp filesep 'TransformParameters.1.txt'];
                At_nr = textread(Trafo_nonrig_result,'%s','delimiter','\n','bufsize',2^18);
                
            end
            
            for i=1:length(for_trafo.trafo_names)
                
                [Fol,dummy] = fileparts(for_trafo.trafo_names{i});
                
                if par.R2D.type~=3
                    FinTrafoN{i} = [Fol filesep 'Final_Trafo.txt'];
                    At{4} = ['(InitialTransformParametersFileName "' for_trafo.trafo_names{i} '")'];
                    suc = EDTI_Library.E_DTI_SMECEPI_write_tra_file(At,FinTrafoN{i});
                    if suc==0
                        EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                        return;
                    end
                else
                    FinTrafoN_R{i} = [Fol filesep 'Final_Trafo_Rigid_step.txt'];
                    FinTrafoN{i} = [Fol filesep 'Final_Trafo.txt'];
                    At{4} = ['(InitialTransformParametersFileName "' for_trafo.trafo_names{i} '")'];
                    suc = EDTI_Library.E_DTI_SMECEPI_write_tra_file(At,FinTrafoN_R{i});
                    if suc==0
                        EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                        return;
                    end
                    At_nr{4} = ['(InitialTransformParametersFileName "' FinTrafoN_R{i} '")'];
                    suc = EDTI_Library.E_DTI_SMECEPI_write_tra_file(At_nr,FinTrafoN{i});
                    if suc==0
                        EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                        return;
                    end
                end
                
            end
            
            
            load(LoadF,'b','bval','info','g','NrB0')
            
            [b, g] = EDTI_Library.E_DTI_reorient_grad_and_b_matrix_rigid_rotation(b,g,Rig);
            
            if isnan(bval)
                diff_model=2;
            else
                diff_model=1;
            end
            
            [dummy, g] =  EDTI_Library.E_DTI_Get_Gradients_and_Bvalue(b, NrB0, diff_model);
            
            parfor i=1:length(FinTrafoN)
                
                EDTI_Library.E_DTI_do_the_SMEC_trafo_step(for_trafo.f_moving{i},for_trafo.dir_temp_i{i},FinTrafoN{i},par)
                
            end
            
            result_names = for_trafo.result_names;
            
            parfor i=1:length(FinTrafoN)
                DWI{i} = EDTI_Library.E_DTI_read_nifti_file(result_names{i});
            end
            
            for i=1:length(DWI)
                if isempty(DWI{i})
                    disp('Errors encountered...')
                    EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                    suc = 0;
                    return;
                end
            end
            
            for i=1:length(DWI)
                DWI{i} = single(DWI{i});
                DWI{i}(DWI{i}==-1000)=nan;
            end
            
            
            dummy = false(size(DWI{1}));
            for i=1:length(DWI)
                dummy = or(dummy,isnan(DWI{i}));
                DWI{i}(isnan(DWI{i}))=0;
                DWI{i}(DWI{i}<0)=0;
            end
            
            if isempty(par.cust_mask.TS)
                mask = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,par.mask_P.TS.NDWI,par.mask_P.TS.DWI,par.mask_P.TS.mfs);
            else
                fn_cm = [FOLD filesep FFO par.cust_mask.TS];
                [mask, VDims, suc] = EDTI_Library.E_DTI_read_nifti_file(fn_cm);
                if suc==0
                    EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                    return;
                end
                mask = mask>0;
            end
            
            mask(dummy)=0;
            
            par_temp = par;
            par_temp.TE = par.TE.TS;
            
            if diff_model==1
                [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par_temp,NrB0);
            elseif diff_model==2
                [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par_temp,NrB0,VDims);
            end
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
            
            g = g(NrB0+1:end,:);
            
            if ~isreal(FA)
                FA = real(FA);
            end
            if ~isreal(FEFA)
                FEFA = real(FEFA);
            end
            if ~isreal(FE)
                FE = real(FE);
            end
            if ~isreal(SE)
                SE = real(SE);
            end
            if ~isreal(eigval)
                eigval = real(eigval);
            end
            MDims = size(mask);
            VDims = VDims_E;
            
            f_out = [par.out_folder filesep FFO par.suff.TS];
            par.Rig = Rig;
            
            max_DWI = 0;
            for i=1:length(DWI)
                max_DWI = max(max_DWI,max(DWI{i}(:)));
            end
            
            if max_DWI<=intmax('int16')
                for i=1:length(DWI)
                    DWI{i} = round(DWI{i});
                    DWI{i} = int16(DWI{i});
                end
            elseif max_DWI<=intmax('uint16')
                for i=1:length(DWI)
                    DWI{i} = round(DWI{i});
                    DWI{i}(DWI{i}<0)=0;
                    DWI{i} = uint16(DWI{i});
                end
            end
            
            
            try
                save(f_out,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
                    'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','par')
            catch me
                EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                disp(me.message)
                return;
            end
            
            if diff_model==2
                save(f_out,'KT','-append')
            end
            
            EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function suc = E_DTI_SMECEPI_make_unity_transf(fn,VDims,MDims)
            
            suc = 1;
            
            VDims = VDims([2 1 3]);
            VDims(3) = -VDims(3);
            
            MDims = MDims([2 1 3]);
            
            ORI = (VDims.*MDims)/2;
            
            c = 0;
            
            c=c+1;L{c} = '(Transform "AffineDTITransform")';
            c=c+1;L{c} = '(NumberOfParameters 12)';
            c=c+1;L{c} = '(TransformParameters 0 0 0 0 0 0 1 1 1 0 0 0)';
            c=c+1;L{c} = '(InitialTransformParametersFileName "NoInitialTransform")';
            c=c+1;L{c} = '(HowToCombineTransforms "Compose")';
            
            c=c+1;L{c} = '// Image specific';
            c=c+1;L{c} = '(FixedImageDimension 3)';
            c=c+1;L{c} = '(MovingImageDimension 3)';
            c=c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c=c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c=c+1;L{c} = ['(Size ' num2str(MDims) ')'];
            c=c+1;L{c} = '(Index 0 0 0)';
            c=c+1;L{c} = ['(Spacing ' num2str(abs(VDims)) ')'];
            c=c+1;L{c} = ['(Origin ' num2str(ORI) ')'];
            c=c+1;L{c} = '(Direction -1 0 0 0 -1 0 0 0 1)';
            c=c+1;L{c} = '(UseDirectionCosines "true")';
            
            c=c+1;L{c} = '// AffineDTITransform specific';
            c=c+1;L{c} = ['(CenterOfRotationPoint ' num2str(VDims/2) ')'];
            c=c+1;L{c} = '(MatrixTranslation 1 0 0 0 1 0 0 0 1 0 0 0)';
            
            c=c+1;L{c} = '// ResampleInterpolator specific';
            c=c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c=c+1;L{c} = '(FinalBSplineInterpolationOrder 1)';
            
            c=c+1;L{c} = '// Resampler specific';
            c=c+1;L{c} = '(Resampler "DefaultResampler")';
            c=c+1;L{c} = '(DefaultPixelValue -1000.000000)';
            c=c+1;L{c} = '(ResultImageFormat "nii")';
            c=c+1;L{c} = '(ResultImagePixelType "float")';
            c=c+1;L{c} = '(CompressResultImage "false")';
            
            
            [fid,mess] = fopen(fn,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' fn])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(fn,'+w -a -s');
            else
                [su, mess] = fileattrib(fn,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' fn])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Creates an affine transformation file for Elastix, tweaked for dMRI data
        function suc = E_DTI_Par_MDC_aff_DTI(par)
            
            suc = 1;
            
            c = 0;
            c = c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c = c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c = c+1;L{c} = '(FixedImageDimension 3)';
            c = c+1;L{c} = '(MovingImageDimension 3)';
            c = c+1;L{c} = '(UseDirectionCosines "true")';
            c = c+1;L{c} = '(Registration "MultiResolutionRegistration")';
            c = c+1;L{c} = '(Interpolator "BSplineInterpolator")';
            c = c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c = c+1;L{c} = '(Resampler "DefaultResampler")';
            c = c+1;L{c} = '(FixedImagePyramid "FixedRecursiveImagePyramid")';
            c = c+1;L{c} = '(MovingImagePyramid "MovingRecursiveImagePyramid")';
            c = c+1;L{c} = '(Optimizer "AdaptiveStochasticGradientDescent")';
            c = c+1;L{c} = '(Transform "AffineDTITransform")';
            c = c+1;L{c} = '(Metric "AdvancedMattesMutualInformation")';
            c = c+1;L{c} = '(AutomaticScalesEstimation "true")';
            c = c+1;L{c} = '(AutomaticTransformInitialization "true")';
            c = c+1;L{c} = '(HowToCombineTransforms "Compose")';
            c = c+1;L{c} = ['(NumberOfHistogramBins ' num2str(par.Hist_bin) ')'];
            c = c+1;L{c} = '(ErodeMask "false")';
            if par.DOF==2
                scales = [-1 -1 -1 ones(1,6)*3.0e+038 -1 -1 -1];
            else
                if par.MDC_constr_fac==0
                    scales = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
                else
                    scales = [-1 -1 -1 10^par.MDC_constr_fac 10^par.MDC_constr_fac ...
                        10^par.MDC_constr_fac 10^par.MDC_constr_fac 10^par.MDC_constr_fac ...
                        10^par.MDC_constr_fac -1 -1 -1];
                end
            end
            c = c+1;L{c} = ['(Scales ' num2str(scales,'%-e ') ')'];
            c = c+1;L{c} = ['(NumberOfResolutions ' num2str(par.Num_Resol) ')'];
            c = c+1;L{c} = ['(MaximumNumberOfIterations ' num2str(par.Num_iter) ')'];
            c = c+1;L{c} = ['(NumberOfSpatialSamples ' num2str(par.Num_samp) ')'];
            c = c+1;L{c} = '(CheckNumberOfSamples "false")';
            c = c+1;L{c} = '(NewSamplesEveryIteration "true")';
            c = c+1;L{c} = '(ImageSampler "RandomCoordinate")';
            if par.Interpol==2
                INTERPOL = 3;
            elseif par.Interpol==1
                INTERPOL = 1;
            end
            % c = c+1;L{c} = ['(BSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(BSplineInterpolationOrder 3)';
            c = c+1;L{c} = ['(FinalBSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(DefaultPixelValue -1000)';
            c = c+1;L{c} = '(WriteResultImage "false")';
            c = c+1;L{c} = '(ResultImagePixelType "float")';
            c = c+1;L{c} = '(ResultImageFormat "nii")';
            c = c+1;L{c} = '(MaximumNumberOfSamplingAttempts 100)';
            
            [fid,mess] = fopen(par.Par_SMEC,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' par.Par_SMEC])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(par.Par_SMEC,'+w -a -s');
            else
                [su, mess] = fileattrib(par.Par_SMEC,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' par.Par_SMEC])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Creates a rigid transformation file for Elastix for EPI purposes
        % Sets Normalized correlation as metric
        function suc = E_DTI_Par_Trafo_DTI_rigid_NC(par)
            
            suc = 1;
            
            c = 0;
            c = c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c = c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c = c+1;L{c} = '(FixedImageDimension 3)';
            c = c+1;L{c} = '(MovingImageDimension 3)';
            c = c+1;L{c} = '(UseDirectionCosines "true")';
            c = c+1;L{c} = '(Registration "MultiResolutionRegistration")';
            c = c+1;L{c} = '(Interpolator "BSplineInterpolator")';
            c = c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c = c+1;L{c} = '(Resampler "DefaultResampler")';
            c = c+1;L{c} = '(FixedImagePyramid "FixedRecursiveImagePyramid")';
            c = c+1;L{c} = '(MovingImagePyramid "MovingRecursiveImagePyramid")';
            c = c+1;L{c} = '(Optimizer "AdaptiveStochasticGradientDescent")';
            c = c+1;L{c} = '(Transform "AffineDTITransform")';
            c = c+1;L{c} = '(Metric "AdvancedNormalizedCorrelation")';
            c = c+1;L{c} = '(AutomaticScalesEstimation "true")';
            c = c+1;L{c} = '(AutomaticTransformInitialization "true")';
            if(isfield(par,'failsafe') && par.failsafe == 1)
                par.Num_Resol = 1;
            else
                c = c+1;L{c} = '(AutomaticTransformInitializationMethod "CenterOfGravity")'; % NEW
            end
            c = c+1;L{c} = '(HowToCombineTransforms "Compose")';
            c = c+1;L{c} = '(ErodeMask "false")';
            scales = [-1 -1 -1 ones(1,6)*3.0e+038 -1 -1 -1];
            c = c+1;L{c} = ['(Scales ' num2str(scales,'%-e ') ')'];
            c = c+1;L{c} = ['(NumberOfHistogramBins ' num2str(par.Hist_bin) ')'];
            % c = c+1;L{c} = '(NumberOfHistogramBins 64)';
            c = c+1;L{c} = ['(NumberOfResolutions ' num2str(par.Num_Resol) ')'];
            % c = c+1;L{c} = '(NumberOfResolutions 2)';
            c = c+1;L{c} = ['(MaximumNumberOfIterations ' num2str(par.Num_iter) ')'];
            % c = c+1;L{c} = '(MaximumNumberOfIterations 1000)';
            c = c+1;L{c} = ['(NumberOfSpatialSamples ' num2str(par.Num_samp) ')'];
            % c = c+1;L{c} = '(NumberOfSpatialSamples 10000)';
            c = c+1;L{c} = '(CheckNumberOfSamples "false")';
            c = c+1;L{c} = '(NewSamplesEveryIteration "true")';
            c = c+1;L{c} = '(ImageSampler "RandomCoordinate")';
            if par.Interpol==2
                INTERPOL = 3;
            elseif par.Interpol==1
                INTERPOL = 1;
            end
            % c = c+1;L{c} = ['(BSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(BSplineInterpolationOrder 3)';
            c = c+1;L{c} = ['(FinalBSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            % c = c+1;L{c} = '(FinalBSplineInterpolationOrder 3)';
            c = c+1;L{c} = '(DefaultPixelValue -1000)';
            c = c+1;L{c} = '(WriteResultImage "false")';
            c = c+1;L{c} = '(ResultImagePixelType "float")';
            c = c+1;L{c} = '(ResultImageFormat "nii")';
            c = c+1;L{c} = '(MaximumNumberOfSamplingAttempts 100)';
            
            [fid,mess] = fopen(par.Par_FN,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' par.Par_FN])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(par.Par_FN,'+w -a -s');
            else
                [su, mess] = fileattrib(par.Par_FN,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' par.Par_FN])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Creates a rigid transformation file for Elastix for EPI purposes
        % Sets mutual information as metric
        function suc = E_DTI_Par_Trafo_DTI_rigid(par)
            
            suc = 1;
            
            c = 0;
            c = c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c = c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c = c+1;L{c} = '(FixedImageDimension 3)';
            c = c+1;L{c} = '(MovingImageDimension 3)';
            c = c+1;L{c} = '(UseDirectionCosines "true")';
            c = c+1;L{c} = '(Registration "MultiResolutionRegistration")';
            c = c+1;L{c} = '(Interpolator "BSplineInterpolator")';
            c = c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c = c+1;L{c} = '(Resampler "DefaultResampler")';
            c = c+1;L{c} = '(FixedImagePyramid "FixedRecursiveImagePyramid")';
            c = c+1;L{c} = '(MovingImagePyramid "MovingRecursiveImagePyramid")';
            c = c+1;L{c} = '(Optimizer "AdaptiveStochasticGradientDescent")';
            c = c+1;L{c} = '(Transform "AffineDTITransform")';
            c = c+1;L{c} = '(Metric "AdvancedMattesMutualInformation")';
            c = c+1;L{c} = '(AutomaticScalesEstimation "true")';
            c = c+1;L{c} = '(AutomaticTransformInitialization "true")';
            if(isfield(par,'failsafe') && par.failsafe == 1)
                par.Num_Resol = 1;
            else
                c = c+1;L{c} = '(AutomaticTransformInitializationMethod "CenterOfGravity")'; % NEW
            end
            c = c+1;L{c} = '(HowToCombineTransforms "Compose")';
            c = c+1;L{c} = '(ErodeMask "false")';
            scales = [-1 -1 -1 ones(1,6)*3.0e+038 -1 -1 -1];
            c = c+1;L{c} = ['(Scales ' num2str(scales,'%-e ') ')'];
            c = c+1;L{c} = ['(NumberOfHistogramBins ' num2str(par.Hist_bin) ')'];
            % c = c+1;L{c} = '(NumberOfHistogramBins 64)';
            c = c+1;L{c} = ['(NumberOfResolutions ' num2str(par.Num_Resol) ')'];
            % c = c+1;L{c} = '(NumberOfResolutions 2)';
            c = c+1;L{c} = ['(MaximumNumberOfIterations ' num2str(par.Num_iter) ')'];
            % c = c+1;L{c} = '(MaximumNumberOfIterations 1000)';
            c = c+1;L{c} = ['(NumberOfSpatialSamples ' num2str(par.Num_samp) ')'];
            % c = c+1;L{c} = '(NumberOfSpatialSamples 10000)';
            c = c+1;L{c} = '(CheckNumberOfSamples "false")';
            c = c+1;L{c} = '(NewSamplesEveryIteration "true")';
            c = c+1;L{c} = '(ImageSampler "RandomCoordinate")';
            if par.Interpol==2
                INTERPOL = 3;
            elseif par.Interpol==1
                INTERPOL = 1;
            end
            % c = c+1;L{c} = ['(BSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(BSplineInterpolationOrder 3)';
            c = c+1;L{c} = ['(FinalBSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            % c = c+1;L{c} = '(FinalBSplineInterpolationOrder 3)';
            c = c+1;L{c} = '(DefaultPixelValue -1000)';
            c = c+1;L{c} = '(WriteResultImage "false")';
            c = c+1;L{c} = '(ResultImagePixelType "float")';
            c = c+1;L{c} = '(ResultImageFormat "nii")';
            c = c+1;L{c} = '(MaximumNumberOfSamplingAttempts 100)';
            
            [fid,mess] = fopen(par.Par_FN,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' par.Par_FN])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(par.Par_FN,'+w -a -s');
            else
                [su, mess] = fileattrib(par.Par_FN,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' par.Par_FN])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Creates a b-spline transformation file for Elastix for EPI purposes
        % Sets Normalized correlation as metric
        function suc = E_DTI_Par_Trafo_DTI_non_rigid_NC(par)
            
            suc = 1;
            
            c = 0;
            c = c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c = c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c = c+1;L{c} = '(FixedImageDimension 3)';
            c = c+1;L{c} = '(MovingImageDimension 3)';
            c = c+1;L{c} = '(UseDirectionCosines "true")';
            c = c+1;L{c} = '(Registration "MultiResolutionRegistration")';
            c = c+1;L{c} = '(Interpolator "BSplineInterpolator")';
            c = c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c = c+1;L{c} = '(Resampler "DefaultResampler")';
            c = c+1;L{c} = '(FixedImagePyramid "FixedRecursiveImagePyramid")';
            c = c+1;L{c} = '(MovingImagePyramid "MovingRecursiveImagePyramid")';
            c = c+1;L{c} = '(Optimizer "AdaptiveStochasticGradientDescent")';
            c = c+1;L{c} = '(Transform "BSplineTransform")';
            c = c+1;L{c} = '(Metric "AdvancedNormalizedCorrelation")';
            c = c+1;L{c} = ['(FinalGridSpacingInPhysicalUnits ' num2str(par.Grid_Spacing) ')'];
            % c = c+1;L{c} = '(FinalGridSpacingInPhysicalUnits 40 40 40)';
            c = c+1;L{c} = '(AutomaticScalesEstimation "true")';
            c = c+1;L{c} = '(AutomaticTransformInitialization "true")';
            c = c+1;L{c} = ['(MovingImageDerivativeScales ' num2str(par.Deriv_Scales) ')'];
            c = c+1;L{c} = '(HowToCombineTransforms "Compose")';
            c = c+1;L{c} = ['(NumberOfHistogramBins ' num2str(par.Hist_bin) ')'];
            % c = c+1;L{c} = '(NumberOfHistogramBins 64)';
            c = c+1;L{c} = '(ErodeMask "false")';
            % c = c+1;L{c} = ['(NumberOfResolutions ' num2str(par.Num_Resol) ')'];
            c = c+1;L{c} = '(NumberOfResolutions 1)';
            c = c+1;L{c} = ['(MaximumNumberOfIterations ' num2str(par.Num_iter) ')'];
            c = c+1;L{c} = ['(NumberOfSpatialSamples ' num2str(par.Num_samp) ')'];
            % c = c+1;L{c} = '(NumberOfResolutions 1)';
            % c = c+1;L{c} = '(MaximumNumberOfIterations 1000)';
            % c = c+1;L{c} = '(NumberOfSpatialSamples 10000)';
            c = c+1;L{c} = '(CheckNumberOfSamples "false")';
            c = c+1;L{c} = '(NewSamplesEveryIteration "true")';
            c = c+1;L{c} = '(ImageSampler "RandomCoordinate")';
            if par.Interpol==2
                INTERPOL = 3;
            elseif par.Interpol==1
                INTERPOL = 1;
            end
            % c = c+1;L{c} = ['(BSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(BSplineInterpolationOrder 3)';
            c = c+1;L{c} = ['(FinalBSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            % c = c+1;L{c} = '(FinalBSplineInterpolationOrder 1)';
            c = c+1;L{c} = '(DefaultPixelValue -1000)';
            c = c+1;L{c} = '(WriteResultImage "false")';
            c = c+1;L{c} = '(ResultImagePixelType "float")';
            c = c+1;L{c} = '(ResultImageFormat "nii")';
            c = c+1;L{c} = '(MaximumNumberOfSamplingAttempts 5)';
            
            [fid,mess] = fopen(par.Par_FN,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' par.Par_FN])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(par.Par_FN,'+w -a -s');
            else
                [su, mess] = fileattrib(par.Par_FN,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' par.Par_FN])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Creates a b-spline transformation file for Elastix for EPI purposes
        % Sets mutual information as metric
        function suc = E_DTI_Par_Trafo_DTI_non_rigid(par)
            
            suc = 1;
            
            c = 0;
            c = c+1;L{c} = '(FixedInternalImagePixelType "float")';
            c = c+1;L{c} = '(MovingInternalImagePixelType "float")';
            c = c+1;L{c} = '(FixedImageDimension 3)';
            c = c+1;L{c} = '(MovingImageDimension 3)';
            c = c+1;L{c} = '(UseDirectionCosines "true")';
            c = c+1;L{c} = '(Registration "MultiResolutionRegistration")';
            c = c+1;L{c} = '(Interpolator "BSplineInterpolator")';
            c = c+1;L{c} = '(ResampleInterpolator "FinalBSplineInterpolator")';
            c = c+1;L{c} = '(Resampler "DefaultResampler")';
            c = c+1;L{c} = '(FixedImagePyramid "FixedRecursiveImagePyramid")';
            c = c+1;L{c} = '(MovingImagePyramid "MovingRecursiveImagePyramid")';
            c = c+1;L{c} = '(Optimizer "AdaptiveStochasticGradientDescent")';
            c = c+1;L{c} = '(Transform "BSplineTransform")';
            c = c+1;L{c} = '(Metric "AdvancedMattesMutualInformation")';
            c = c+1;L{c} = ['(FinalGridSpacingInPhysicalUnits ' num2str(par.Grid_Spacing) ')'];
            % c = c+1;L{c} = '(FinalGridSpacingInPhysicalUnits 40 40 40)';
            c = c+1;L{c} = '(AutomaticScalesEstimation "true")';
            c = c+1;L{c} = '(AutomaticTransformInitialization "true")';
            c = c+1;L{c} = ['(MovingImageDerivativeScales ' num2str(par.Deriv_Scales) ')'];
            c = c+1;L{c} = '(HowToCombineTransforms "Compose")';
            c = c+1;L{c} = ['(NumberOfHistogramBins ' num2str(par.Hist_bin) ')'];
            % c = c+1;L{c} = '(NumberOfHistogramBins 64)';
            c = c+1;L{c} = '(ErodeMask "false")';
            % c = c+1;L{c} = ['(NumberOfResolutions ' num2str(par.Num_Resol) ')'];
            c = c+1;L{c} = '(NumberOfResolutions 1)';
            c = c+1;L{c} = ['(MaximumNumberOfIterations ' num2str(par.Num_iter) ')'];
            c = c+1;L{c} = ['(NumberOfSpatialSamples ' num2str(par.Num_samp) ')'];
            % c = c+1;L{c} = '(NumberOfResolutions 1)';
            % c = c+1;L{c} = '(MaximumNumberOfIterations 1000)';
            % c = c+1;L{c} = '(NumberOfSpatialSamples 10000)';
            c = c+1;L{c} = '(CheckNumberOfSamples "false")';
            c = c+1;L{c} = '(NewSamplesEveryIteration "true")';
            c = c+1;L{c} = '(ImageSampler "RandomCoordinate")';
            if par.Interpol==2
                INTERPOL = 3;
            elseif par.Interpol==1
                INTERPOL = 1;
            end
            % c = c+1;L{c} = ['(BSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            c = c+1;L{c} = '(BSplineInterpolationOrder 3)';
            c = c+1;L{c} = ['(FinalBSplineInterpolationOrder ' num2str(INTERPOL) ')'];
            % c = c+1;L{c} = '(FinalBSplineInterpolationOrder 1)';
            c = c+1;L{c} = '(DefaultPixelValue -1000)';
            c = c+1;L{c} = '(WriteResultImage "false")';
            c = c+1;L{c} = '(ResultImagePixelType "float")';
            c = c+1;L{c} = '(ResultImageFormat "nii")';
            c = c+1;L{c} = '(MaximumNumberOfSamplingAttempts 5)';
            
            [fid,mess] = fopen(par.Par_FN,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' par.Par_FN])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(par.Par_FN,'+w -a -s');
            else
                [su, mess] = fileattrib(par.Par_FN,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' par.Par_FN])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
            
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function E_DTI_output_TR_files(DWI,Vox_S,f_moving,par)
            
            if par.DOF~=1
                if par.Regul==1
                    DWI = EDTI_Library.E_DTI_Anisotr_Gauss_smoothing(DWI,Vox_S,2*max(Vox_S),[3 3 3]);
                end
            end
            
            EDTI_Library.E_DTI_write_nifti_file(DWI,Vox_S,f_moving);
        end
        
        % From ExploreDTI: perform registration - using Elastix
        function E_DTI_do_the_SMEC_reg_step(f_moving,dir_temp_i,par)
            
            if ispc
                
                sys_c = ['"' par.E_path '"' ' -f ' '"' par.fixed_SMEC_fn '"'...
                    ' -fMask ' '"' par.fixed_SMEC_fn_mask '"' ...
                    ' -m ' '"' f_moving '"'...
                    ' -out ' '"' dir_temp_i '"' ...
                    ' -p ' '"' par.Par_SMEC '"'];
                [r,st]=system(sys_c);
                
            else
                
                %      par.fixed_SMEC_fn = replace(par.fixed_SMEC_fn,' ','\ ');
                %      par.fixed_SMEC_fn_mask = replace(par.fixed_SMEC_fn_mask,' ','\ ');
                %      f_moving = replace(f_moving,' ','\ ');
                %      dir_temp_i = replace(dir_temp_i,' ','\ ');
                %      par.Par_SMEC = replace(par.Par_SMEC,' ','\ ');
                
                [dir_e,el] = fileparts(par.E_path);
                %     sys_c = ['./' el ' -f ' par.fixed_SMEC_fn ...
                %         ' -fMask ' par.fixed_SMEC_fn_mask ...
                %         ' -m ' f_moving ...
                %         ' -out ' dir_temp_i ...
                %         ' -p ' par.Par_SMEC];
                if(ismac)
                    sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' par.fixed_SMEC_fn ...
                        ' -fMask ' par.fixed_SMEC_fn_mask ...
                        ' -m ' f_moving ...
                        ' -out ' dir_temp_i ...
                        ' -p ' par.Par_SMEC];
                else
                    sys_c = ['LD_LIBRARY_PATH="' dir_e '" "' par.E_path '" -f ' par.fixed_SMEC_fn ...
                        ' -fMask ' par.fixed_SMEC_fn_mask ...
                        ' -m ' f_moving ...
                        ' -out ' dir_temp_i ...
                        ' -p ' par.Par_SMEC];
                end
                %     tedi = cd(dir_e);
                [r,st]=system(sys_c);
                %     cd(tedi)
            end
            
            if r~=0
                disp('Error(s) during motion/distortion correction:')
                disp(' ')
                disp(st);
                disp(' ')
                disp('For a solution, check the forum...')
            end
            
        end
        
        % From ExploreDTI: perform co-registration - using Elastix/Transformix
        function E_DTI_do_the_SMEC_trafo_step(f_moving,dir_temp_i,tra_file,par)
            
            if ispc
                
                sys_c = ['"' par.T_path '"' ' -in ' '"' f_moving '"'...
                    ' -out ' '"' dir_temp_i '"' ...
                    ' -tp ' '"' tra_file '"'];
                [r,st]=system(sys_c);
                
            else
                
                %     dir_temp_i = replace(dir_temp_i,' ','\ ');
                %     dir_temp_i = ['""' dir_temp_i '""'];
                %     tra_file = replace(tra_file,' ','\ ');
                %     f_moving = replace(f_moving,' ','\ ');
                
                [dir_e,tra] = fileparts(par.T_path);
                %     sys_c = ['./' tra ' -in ' f_moving ...
                %         ' -out ' dir_temp_i ...
                %         ' -tp ' tra_file];
                if(ismac)
                    sys_c = ['DYLD_LIBRARY_PATH="' dir_e '" "' par.T_path '" -in ' f_moving ...
                        ' -out ' dir_temp_i ...
                        ' -tp ' tra_file];
                else
                    sys_c = ['LD_LIBRARY_PATH="' dir_e '" "' par.T_path '" -in ' f_moving ...
                        ' -out ' dir_temp_i ...
                        ' -tp ' tra_file];
                end
                %     tedi = cd(dir_e);
                [r,st]=system(sys_c);
                %     cd(tedi)
                
            end
            
            if r~=0
                disp('Error(s) during motion/distortion correction:')
                disp(' ')
                disp(st);
                disp(' ')
                disp('For a solution, check the forum...')
            end
            
        end
        
        % From ExploreDTI: read transformation matrices to perform b-matrix
        % rotation
        function Trafo = E_DTI_Tra_Par_2_Tra_Mat(Par)
            
            % Example
            
            % Par = ...
            %    [10.1776  -10.3675    8.4619;
            %     2.5370   -6.4122    5.1677;
            %     0.9089    1.1103    0.9521;
            %    -0.0108    0.0197   -0.0304];
            
            % with the following definition:
            
            % Par =    [RotX   RotY   RotZ   ;
            %           Tx     Ty     Tz     ;
            %           ScaleX ScaleY ScaleZ ;
            %           SkewX  SkewY  SkewZ] ;
            
            % results in:
            
            % Trafo = [0.8869    0.1550    0.1613    2.5370;
            %         -0.1429    1.0698    0.1673   -6.4122;
            %         -0.1392   -0.2557    0.9233    5.1677];
            
            
            tx = Par(2,1);
            ty = Par(2,2);
            tz = Par(2,3);
            
            rx = deg2rad(Par(1,1));
            ry = deg2rad(Par(1,2));
            rz = deg2rad(Par(1,3));
            
            sx = Par(3,1);
            sy = Par(3,2);
            sz = Par(3,3);
            
            gx = Par(4,1);
            gy = Par(4,2);
            gz = Par(4,3);
            
            T = [1 0 0 tx;
                0 1 0 ty;
                0 0 1 tz;
                0 0 0 1];
            
            Rx = [1 0 0 0;
                0 cos(rx) sin(rx) 0;
                0 -sin(rx) cos(rx) 0;
                0 0 0 1];
            
            Ry = [cos(ry) 0 -sin(ry) 0;
                0 1 0 0;
                sin(ry) 0 cos(ry) 0;
                0 0 0 1];
            
            Rz = [cos(rz) sin(rz) 0 0;
                -sin(rz) cos(rz) 0 0;
                0 0 1 0;
                0 0 0 1];
            
            R = Rx*Ry*Rz;
            
            Gx = [1 0 gx 0;
                0 1 0 0;
                0 0 1 0;
                0 0 0 1];
            
            Gy = [1 0 0 0;
                gy 1 0 0;
                0 0 1 0;
                0 0 0 1];
            
            Gz = [1 0 0 0;
                0 1 0 0;
                0 gz 1 0;
                0 0 0 1];
            
            G = Gx*Gy*Gz;
            
            S = [sx 0 0 0;
                0 sy 0 0;
                0 0 sz 0;
                0 0 0 1];
            
            
            Trafo = T*R*G*S;
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function suc = E_DTI_SMECEPI_write_tra_file(L,FN)
            
            suc = 1;
            
            [fid,mess] = fopen(FN,'w+t');
            
            if fid==-1
                suc=0;
                disp(['Problem with making file: ' FN])
                disp(['Error message: ' mess])
                return;
            end
            
            for i=1:length(L)
                fprintf(fid, '%s\n', L{i});
            end
            
            try
                fclose(fid);
            catch me
                suc = 0;
                disp(me.message)
                return;
            end
            
            if ispc
                [su, mess] = fileattrib(FN,'+w -a -s');
            else
                [su, mess] = fileattrib(FN,'+w +x','a');
            end
            
            if su==0
                disp(['Problem with file: ' FN])
                disp(['Error message: ' mess])
                suc = 0;
                return;
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Concatenates rotation matrices for b-matrix rotation
        function Trafo = E_DTI_Convert_aff_par_2_aff_transf_matrix(Par)
            
            % Example
            
            % Par = ...
            %    [10.1776  -10.3675    8.4619;
            %     2.5370   -6.4122    5.1677;
            %     0.9089    1.1103    0.9521;
            %    -0.0108    0.0197   -0.0304];
            
            % with the following definition:
            
            % Par =    [RotX   RotY   RotZ   ;
            %           Tx     Ty     Tz     ;
            %           ScaleX ScaleY ScaleZ ;
            %           SkewX  SkewY  SkewZ] ;
            
            % results in:
            
            % Trafo = [0.8869    0.1550    0.1613    2.5370;
            %         -0.1429    1.0698    0.1673   -6.4122;
            %         -0.1392   -0.2557    0.9233    5.1677];
            
            
            tx = Par(2,1);
            ty = Par(2,2);
            tz = Par(2,3);
            
            rx = deg2rad(Par(1,1));
            ry = deg2rad(Par(1,2));
            rz = deg2rad(Par(1,3));
            
            sx = Par(3,1);
            sy = Par(3,2);
            sz = Par(3,3);
            
            gx = Par(4,1);
            gy = Par(4,2);
            gz = Par(4,3);
            
            T = [1 0 0 tx;
                0 1 0 ty;
                0 0 1 tz;
                0 0 0 1];
            
            Rx = [1 0 0 0;
                0 cos(rx) sin(rx) 0;
                0 -sin(rx) cos(rx) 0;
                0 0 0 1];
            
            Ry = [cos(ry) 0 -sin(ry) 0;
                0 1 0 0;
                sin(ry) 0 cos(ry) 0;
                0 0 0 1];
            
            Rz = [cos(rz) sin(rz) 0 0;
                -sin(rz) cos(rz) 0 0;
                0 0 1 0;
                0 0 0 1];
            
            R = Rx*Ry*Rz;
            
            Gx = [1 0 gx 0;
                0 1 0 0;
                0 0 1 0;
                0 0 0 1];
            
            Gy = [1 0 0 0;
                gy 1 0 0;
                0 0 1 0;
                0 0 0 1];
            
            Gz = [1 0 0 0;
                0 1 0 0;
                0 gz 1 0;
                0 0 0 1];
            
            G = Gx*Gy*Gz;
            
            S = [sx 0 0 0;
                0 sy 0 0;
                0 0 sz 0;
                0 0 0 1];
            
            
            Trafo = T*R*G*S;
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Performs the b-matrix rotation. 
        function [b, g] = E_DTI_reorient_grad_and_b_matrix_rigid_rotation(b,g,Rot)
            M = zeros(4,3);
            M(1,:) = [Rot{1} Rot{2} Rot{3}]*(180/pi);
            M(3,:) = 1;
            World_Trafo = EDTI_Library.E_DTI_Convert_aff_par_2_aff_transf_matrix(M);
            R = World_Trafo(1:3,1:3);
        
            b_old = b;
            g_old = g;
            
            for i=1:size(b,1)
                
                B = [b_old(i,1) b_old(i,2)/2 b_old(i,3)/2;...
                    b_old(i,2)/2 b_old(i,4) b_old(i,5)/2;...
                    b_old(i,3)/2 b_old(i,5)/2 b_old(i,6)];
                B_rot = R*B*(R');
                
                b(i,:) = [B_rot(1,1) 2*B_rot(1,2) 2*B_rot(1,3)...
                    B_rot(2,2) 2*B_rot(2,3) B_rot(3,3)];
                
            end
            
            
            for i=1:size(g,1)
                
                g(i,:) = g_old(i,:)*R';
                
            end
        end
        
        % From ExploreDTI: small helper function
        function r = deg2rad(d)
            
            r = (pi*d)/180;
            
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function [suc, FN] = E_DTI_SMECEPI_check_data_stuff(f_in,par)
            
            suc = 1;
            
            if par.R2D.type==0
                FN = par.R2D.FN;
                return;
            elseif par.R2D.type==1
                FN = par.R2D.FN;
            else
                %     [Fol,Fil] = fileparts(f_in);
                [Fol,Fil] = fileparts(par.R2D.FN);
                if(~isempty(Fol))
                    %         FN = [Fol filesep Fil par.R2D.FN];
                    FN = par.R2D.FN;
                else
                    FN = fullfile(pwd,par.R2D.FN);
                end
            end
            
            [D, VDims, suc] = EDTI_Library.E_DTI_read_nifti_file(FN);
            
            if suc==0
                disp('Processing stopped for file:')
                disp(f_in)
                return;
            else
                if ndims(D)~=3
                    disp('Processing stopped for file:')
                    disp(f_in)
                    disp('The data in this file should be 3D!')
                    suc=0;
                    return;
                end
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function E_DTI_remove_temp_f(dir_temp)
            
            % return;
            
            warning off all
            suc = rmdir(dir_temp,'s');
            warning on all
            
            if suc==0
                disp('Could not remove temp folder:')
                disp(dir_temp)
                disp('Remove it manually after processing has finished...')
            end
        end
        
        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        function [bval, grad] =  E_DTI_Get_Gradients_and_Bvalue(bmat, NrB0, model)
            
            bmat = double(bmat);
            if model==1
                bval = mean(sum(bmat(NrB0+1:end,[1 4 6]),2));
            else
                bval = nan;
            end
            
            % x~=0
            BSign_x = sign(sign(bmat(:,1:3)) + 0.0001);   % Adding 0.0001 avoids getting zeros here
            
            % y~=0
            BSign_y = sign(sign(bmat(:,[2 4 5])) + 0.0001);   % Adding 0.0001 avoids getting zeros here
            
            Bo = bmat(:,1)==0;
            
            BSign = BSign_x;
            BSign([Bo Bo Bo]) = BSign_y([Bo Bo Bo]);
            
            grad = BSign .* sqrt(bmat(:,[1 4 6]));
            grad_n = sqrt(sum(grad.*grad,2));
            grad = grad./[grad_n grad_n grad_n];
            
            if ~isreal(grad)
                %     disp('Warning: incorrect values of the B-matrix were encountered!')
                grad = real(grad);
                grad_n = sqrt(sum(grad.*grad,2));
                grad = grad./[grad_n grad_n grad_n];
            end
            
            grad(isnan(grad))=0;
        end
        
        % From ExploreDTI: Anisotropic gaussian smoothing
        function B = E_DTI_Anisotr_Gauss_smoothing(A,VDims,FWHM,siz)
            
            str_cl = class(A);
            
            if FWHM<eps
                B = A;
                return;
            end
            
            sigma_sm = FWHM/(2*sqrt(2*log(2)));
            sigma_av = (3*max(VDims))/(2*sqrt(2*log(2)));
            
            A = double(A);
            min_A = min(A(:));
            max_A = max(A(:));
            
            A = (A-min_A)/(max_A-min_A);
            A(isnan(A)) = 0;
            
            ps = (siz-1)/2;
            ps2 = (siz+1)/2;
            A = padarray(A,ps,'replicate');
            
            [y,x,z] = meshgrid(-ps(1):ps(1),-ps(2):ps(2),-ps(3):ps(3));
            
            D = [x(:)*VDims(1) y(:)*VDims(2) z(:)*VDims(3)];
            
            Frow = convn(A, x*(1/VDims(1)), 'same');
            Fcol = convn(A, y*(1/VDims(2)), 'same');
            Fsli = convn(A, z*(1/VDims(3)), 'same');
            
            S{1} = Frow.*Frow;
            S{2} = Frow.*Fcol;
            S{3} = Frow.*Fsli;
            S{4} = Fcol.*Fcol;
            S{5} = Fcol.*Fsli;
            S{6} = Fsli.*Fsli;
            
            GM = sqrt(S{1} + S{4} + S{6});
            
            for i=1:6
                S{i} = EDTI_Library.E_DTI_Isotropic_Gauss_smoothing(S{i}, siz, sigma_av, VDims);
            end
            
            [FEFA, FA, FE, SE, eigval, Bm] = EDTI_Library.E_DTI_eigensystem_analytic(S);
            TE = cross(FE,SE,4);
            
            N{1} = TE;
            N{2} = SE;
            N{3} = FE;
            
            for i=1:3
                N{i} = N{i}(ps2(1):end+1-ps2(1),...
                    ps2(2):end+1-ps2(2),...
                    ps2(3):end+1-ps2(3),:);
            end
            
            clear TE SE FE;
            
            L1 = eigval(:,:,:,1);
            L2 = eigval(:,:,:,2);
            L3 = eigval(:,:,:,3);
            
            SL = L1+L2+L3;
            SL(SL==0)=1;
            
            a12 = (L2-L3)./SL;
            a13 = (L1-L3)./SL;
            
            C = (1-a12-a13).*GM;
            
            sig{1} = sigma_sm./(1 + C);
            sig{2} = sig{1}.*(1 - 2*a12);
            sig{3} = sig{1}.*(1-a12-a13);
            
            for i=1:3
                sig{i} = sig{i}(ps2(1):end+1-ps2(1),...
                    ps2(2):end+1-ps2(2),...
                    ps2(3):end+1-ps2(3));
            end
            
            L = size(D,1);
            
            Tot_Weight = 0;
            Tot_Sum = 0;
            
            for i=1:L
                
                Dum = A(ps2(1)+x(i):end+1-ps2(1)+x(i),...
                    ps2(2)+y(i):end+1-ps2(2)+y(i),...
                    ps2(3)+z(i):end+1-ps2(3)+z(i));
                
                DumV = 0;
                
                for k=1:3
                    DumV = DumV + ((D(i,1)*N{k}(:,:,:,1) + D(i,2)*N{k}(:,:,:,2) + D(i,3)*N{k}(:,:,:,3)).^2)./((sig{k}).^2);
                end
                
                Weight = exp((-1/2)*DumV);
                
                Tot_Sum = Tot_Sum + Dum.*Weight;
                
                Tot_Weight = Tot_Weight + Weight;
                
            end
            
            B = Tot_Sum./Tot_Weight;
            B = B*(max_A-min_A) + min_A;
            B(isnan(B))=0;
            
            if strcmpi(str_cl,'single')
                B = single(B);
            elseif strcmpi(str_cl,'int16')
                B = round(B);
                B(B>intmax('int16'))=intmax('int16');
                B(B<intmin('int16'))=intmin('int16');
            elseif strcmpi(str_cl,'uint16')
                B = round(B);
                B(B>intmax('uint16'))=intmax('uint16');
                B(B<intmin('uint16'))=intmin('uint16');
            elseif strcmpi(str_cl,'int8')
                B = round(B);
                B(B>intmax('int8'))=intmax('int8');
                B(B<intmin('int8'))=intmin('int8');
            elseif strcmpi(str_cl,'uint8')
                B = round(B);
                B(B>intmax('uint8'))=intmax('uint8');
                B(B<intmin('uint8'))=intmin('uint8');
            end
            
        end
        
        % From ExploreDTI: Isotropic gaussian smoothing
        function D = E_DTI_Isotropic_Gauss_smoothing(data, sz, arg, vox)
            
            if arg==0
                D=data;
                return;
            end
            
            maskn = isnan(data);
            data(maskn) = 0;
            
            if length(sz)==1
                sz = [sz sz sz];
            end
            
            sz = sz(:)';
            padSize = (sz-1)/2;
            smooth = gaussian3(sz,arg,vox);
            
            D=convn(EDTI_Library.padreplicate(data,padSize),smooth, 'valid');
            D(maskn)=nan;
        end
        
        % From ExploreDTI: helper for isotropic gaussian smoothing
        function h = gaussian3(siz, std, vox)
            
            [x,y,z] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2, -(siz(1)-1)/2:(siz(1)-1)/2, -(siz(3)-1)/2:(siz(3)-1)/2);
            x = x*vox(2);
            y = y*vox(1);
            z = z*vox(3);
            h = exp(-(x.*x + y.*y + z.*z)/(2*std*std));
            h = h/sum(h(:));
        end
        
        % From ExploreDTI: helper for isotropic gaussian smoothing
        function b=padreplicate(a, padSize)
            
            numDims = length(padSize);
            idx = cell(numDims,1);
            for k = 1:numDims
                M = size(a,k);
                onesVector = ones(1,padSize(k));
                idx{k} = [onesVector 1:M M*onesVector];
            end
            b = a(idx{:});
        end
        
        % From ExploreDTI: helper function for DKI
        function MK = E_DTI_Mean_Kurtosis_c(KT,DT)
            
            mask = ~isnan(KT{1});
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i})';
            end
            
            for i=[1 4 6]
                %     DT{i}(DT{i}<=0) = 0;
                DT{i}(DT{i}<0) = abs(DT{i}(DT{i}<0));
            end
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            p = which('Grad_dirs_1024.txt');
            g = textread(p);
            
            
            MDsq = repmat(single(nan),size(KT{1}));
            MK = repmat(single(nan),size(KT{1}));
            MK(mask)=0;
            MDsq(mask) = ((DT{1}+DT{4}+DT{6})/3).^2;
            
            dt = [DT{1}; DT{4}; DT{6}; DT{2}; DT{3}; DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i})';
            end
            
            %-
            for i = [1 2 3 10 11 12]
                %     KT{i}(KT{i}<=0) = 0;
                KT{i}(KT{i}<=0) = abs(KT{i}(KT{i}<=0));
            end
            
            kt = [KT{1}; KT{2}; KT{3}; KT{4}; KT{5}; KT{6}; ...
                KT{7}; KT{8}; KT{9}; KT{10}; KT{11}; KT{12}; ...
                KT{13}; KT{14}; KT{15}];
            
            clear KT;
            
            cn = MDsq(mask);
            cn(:) = 0;
            V_MK = cn;
            
            VMDsq = MDsq(mask);
            
            for i=1:size(g,1)
                
                dum_adc = (B(i,:)*dt).^2;
                dum_kt = A(i,:)*kt;
                
                dum_adc = dum_adc(:);
                dum_kt = dum_kt(:);
                
                M_c = dum_kt(:)>0;
                cn(M_c) = cn(M_c)+1;
                
                %     MK(mask(M_c)) = MK(mask) + (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
                V_MK(M_c) = V_MK(M_c) + (VMDsq(M_c).*dum_kt(M_c))./dum_adc(M_c);
            end
            
            V_MK = V_MK./cn;
            
            MK(mask) = V_MK(:);
            
            [dummy, MK] = EDTI_Library.E_DTI_local_outlier_correction_avg(MK,mask);
            
            % y = prctile(MK(mask),[1 99]);
            % MK = EDTI_Library.E_DTI_outlier_smoothing(MK,mask,y);
        end
        
        % From ExploreDTI: helper function for DKI
        function MK = E_DTI_Mean_Kurtosis(KT,DT)
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            global MRIToolkit;
            if(isfield(MRIToolkit,'DKI_cleanup') && MRIToolkit.DKI_cleanup == true)
                MK = EDTI_Library.E_DTI_Mean_Kurtosis_c(KT,DT);
                return;
            end
            
            p = which('Grad_dirs_1024.txt');
            g = textread(p);
            
            mask = ~isnan(KT{1});
            
            MDsq = repmat(single(nan),size(KT{1}));
            MK = repmat(single(nan),size(KT{1}));
            MK(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i})';
            end
            
            dt = [DT{1}; DT{4}; DT{6}; DT{2}; DT{3}; DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i})';
            end
            
            kt = [KT{1}; KT{2}; KT{3}; KT{4}; KT{5}; KT{6}; ...
                KT{7}; KT{8}; KT{9}; KT{10}; KT{11}; KT{12}; ...
                KT{13}; KT{14}; KT{15}];
            
            clear KT;
            
            for i=1:size(g,1)
                
                dum_adc = (B(i,:)*dt).^2;
                dum_kt = A(i,:)*kt;
                
                MK(mask) = MK(mask) + (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
            end
            
            MK(mask) = MK(mask)/size(g,1);
            
            % y = prctile(MK(mask),[1 99]);
            % MK = EDTI_Library.E_DTI_outlier_smoothing(MK,mask,y);
        end
        
        % From ExploreDTI: helper function for DKI
        function [BK, AKA] = E_DTI_local_outlier_correction_avg(AK,mask)
            
            kappa = 1.5;
            
            % h = EDTI_Library.E_DTI_gaussian3([3 3 3], 1, VDims);
            AK = real(AK);
            AKA = AK;
            BK = false(size(AK));
            
            for i=1:size(AK,1)
                for j=1:size(AK,2)
                    for k=1:size(AK,3)
                        
                        if mask(i,j,k)==1
                            D = [];
                            G = AK(i,j,k);
                            
                            for r=-1:1:1
                                for s=-1:1:1
                                    for t=-1:1:1
                                        
                                        u = i+r;
                                        v = j+s;
                                        w = k+t;
                                        
                                        if u>0 && v>0 && w>0 && ...
                                                u<size(AK,1)+1 && ...
                                                v<size(AK,2)+1 && ...
                                                w<size(AK,3)+1
                                            
                                            if mask(u,v,w)==1
                                                
                                                D = [D; AK(u,v,w)];
                                                
                                            end
                                        end
                                    end
                                end
                            end
                            
                            if length(D)>10
                                try
                                    y = prctile(D,[25 75]);
                                catch err
                                    disp(err)
                                end
                                
                                IQR = (y(2)-y(1));
                                
                                if G < y(1)-kappa*IQR || G > y(2)+kappa*IQR
                                    
                                    BK(i,j,k) = 1;
                                    
                                end
                                
                            end
                        end
                        
                    end
                end
            end
            
            
            for i=1:size(AK,1)
                for j=1:size(AK,2)
                    for k=1:size(AK,3)
                        
                        if BK(i,j,k)==1
                            D = [];
                            
                            for r=-1:1:1
                                for s=-1:1:1
                                    for t=-1:1:1
                                        
                                        u = i+r;
                                        v = j+s;
                                        w = k+t;
                                        
                                        if u>0 && v>0 && w>0 && ...
                                                u<size(AK,1)+1 && ...
                                                v<size(AK,2)+1 && ...
                                                w<size(AK,3)+1
                                            
                                            if mask(u,v,w)==1 && BK(u,v,w)~=1
                                                
                                                D = [D; AK(u,v,w)];
                                                
                                            end
                                        end
                                    end
                                end
                            end
                            
                            AKA(i,j,k) = mean(D);
                            
                        end
                        
                    end
                end
            end
        end
        
        % From ExploreDTI: Helper function
        function B = E_DTI_DT_cell2mat(A)
            
            if iscell(A)
                B = repmat(A{1},[1 1 1 9]);
                B(:,:,:,2) = A{2};
                B(:,:,:,3) = A{3};
                B(:,:,:,4) = A{2};
                B(:,:,:,5) = A{4};
                B(:,:,:,6) = A{5};
                B(:,:,:,7) = A{3};
                B(:,:,:,8) = A{5};
                B(:,:,:,9) = A{6};
            else
                B=A;
            end
        end
        
        % From ExploreDTI: read .nii files in the ExploreDTI format
        function [I, VDims, suc, hdr] = E_DTI_read_nifti_file(fn)
            
            if(exist(fn,'file') < 1 && strcmp(fn(end-3:end),'.nii'))
                fn = [fn '.gz'];
            end
            
            % This is a quick-and-dirty *.nii file reader...
            
            [P,N,EXT] = fileparts(fn);
            temp_file = false;
            if strcmpi(EXT,'.gz')
                gzipped = 1;
                while(true)
                    temp_path = fullfile(tempdir,['mrtd_gunzip_' num2str(randi(100000))]);
                    if(exist(temp_path,'dir') < 1)
                        mkdir(temp_path)
                        temp_file = true;
                        break
                    end
                end
                fn = gunzip(fn,temp_path);
                [P2,N2,EXT2] = fileparts(fn{1});
                if isempty(P2)
                    P2 = pwd;
                    fn = [P2 filesep fn{1}];
                else
                    fn = fn{1};
                end
            else
                gzipped = 0;
            end
            
            hdr = load_untouch_header_only(fn);
            
            I = [];
            VDims = [];
            suc = 1;
            
            [fid,mess] = fopen(fn, 'r');
            if fid == -1
                disp('Could not open file:')
                disp(fn)
                disp(['Error message: ' mess])
                suc=0;
                return;
            end
            
            try
                
                header.sizeof_hdr = fread(fid, 1, 'int32');
                header.data_type = fread(fid, 10, 'uint8');
                header.db_name = fread(fid, 18, 'uint8');
                header.extents = fread(fid, 1, 'int32');
                header.session_error = fread(fid, 1, 'int16');
                header.regular = fread(fid, 1, 'uint8');
                header.dim_info = fread(fid, 1, 'uint8');
                header.dim = fread(fid, 8, 'int16');
                header.intent_p1 = fread(fid, 1, 'float32');
                header.intent_p2 = fread(fid, 1, 'float32');
                header.intent_p3 = fread(fid, 1, 'float32');
                header.intent_code = fread(fid, 1, 'int16');
                header.datatype = fread(fid, 1, 'int16');
                header.bitpix = fread(fid, 1, 'int16');
                header.slice_start = fread(fid, 1, 'int16');
                header.pixdim = fread(fid, 8, 'float32');
                header.vox_offset = fread(fid, 1, 'float32');
                header.scl_slope = fread(fid, 1, 'float32');
                header.scl_inter = fread(fid, 1, 'float32');
                header.slice_end = fread(fid, 1, 'int16');
                header.slice_code = fread(fid, 1, 'uint8');
                header.xyzt_units = fread(fid, 1, 'uint8');
                header.cal_max = fread(fid, 1, 'float32');
                header.cal_min = fread(fid, 1, 'float32');
                header.slice_duration = fread(fid, 1, 'float32');
                header.toffset = fread(fid, 1, 'float32');
                header.glmax = fread(fid, 1, 'int32');
                header.glmin = fread(fid, 1, 'int32');
                header.descrip = fread(fid, 80, 'uint8');
                header.aux_file = fread(fid, 24, 'uint8');
                header.qform_code = fread(fid, 1, 'int16');
                header.sform_code = fread(fid, 1, 'int16');
                header.quatern_b = fread(fid, 1, 'float32');
                header.quatern_c = fread(fid, 1, 'float32');
                header.quatern_d = fread(fid, 1, 'float32');
                header.qoffset_x = fread(fid, 1, 'float32');
                header.qoffset_y = fread(fid, 1, 'float32');
                header.qoffset_z = fread(fid, 1, 'float32');
                header.srow_x = fread(fid, 4, 'float32');
                header.srow_y = fread(fid, 4, 'float32');
                header.srow_z = fread(fid, 4, 'float32');
                header.intent_name = fread(fid, 16, 'uint8');
                header.magic = fread(fid, 8, 'char');
                
                siz_hdr = ftell(fid) - 4;
                if (siz_hdr ~= header.sizeof_hdr)
                    fclose(fid);
                    suc=0;
                    return;
                end
                
                fseek(fid,header.vox_offset,-1);
                
                if (header.bitpix == 24)
                    
                    if header.dim(1) > 3
                        suc=0;
                        return;
                    end
                    
                    dimx     = header.dim(2);
                    dimy     = header.dim(3);
                    n_slices = header.dim(4);
                    
                    I = fread(fid, 3*dimx*dimy*n_slices, '*uint8');
                    I = reshape(I, 3, dimx, dimy, n_slices);
                else
                    
                    dimx = header.dim(2);
                    dimy = header.dim(3);
                    n_slices = header.dim(4);
                    n_dir = header.dim(5);
                    
                    switch (header.datatype)
                        case 1280
                            datatype = 'uint64';
                        case 1024
                            datatype = 'int64';
                        case 512
                            datatype = 'uint16';
                        case 256
                            datatype = 'int8';
                        case 64
                            datatype = 'float64';
                        case 16
                            datatype = 'float32';
                        case 8
                            datatype = 'int32';
                        case 4
                            datatype = 'int16';
                        case 2
                            datatype = 'uint8';
                        case 1
                            datatype = 'bit1';
                        otherwise
                            suc=0;
                            return;
                    end
                    
                    I = fread(fid, dimx*dimy*n_slices*n_dir, ['*' datatype]);
                    I = reshape(I, dimx, dimy, n_slices, n_dir);
                end
                
                fclose(fid);
                
                VDims = header.pixdim(2:4)';
                
                if ndims(I)==4
                    I = permute(I,[2 1 3 4]);
                elseif ndims(I)==3 || ndims(I)==2
                    I = permute(I,[2 1 3]);
                elseif ndims(I)==5
                    I = permute(I,[2 1 3 4 5]);
                end
                VDims = VDims([2 1 3]);
                
                I = flipdim(I,1);
                I = flipdim(I,2);
                
                if gzipped==1 || temp_file == true
                    delete(fn);
                    rmdir(temp_path);
                end
                
                
            catch me
                suc = 0;
                disp(me.message)
            end
        end
        
        % From ExploreDTI: write .nii files in the ExploreDTI format
        function E_DTI_write_nifti_file(I, VDims, file_name)
            global MRIToolkit;
            if(isfield(MRIToolkit,'EnforceNiiGz') && MRIToolkit.EnforceNiiGz == true)
                if(strcmp(file_name(end-3:end),'.nii'))
                    file_name = [file_name '.gz'];
                end
            end
            % This is a quick-and-dirty *.nii file writer... Check for left/right
            % flipping, axis permutations, and voxel size first!
            
            [Pat,FN,EXT] = fileparts(file_name);
            if isempty(Pat)
                file_name = [pwd filesep file_name];
            end
            
            [Pat,FN,EXT] = fileparts(file_name);
            
            if strcmpi(EXT,'.gz')
                file_name = [Pat filesep FN];
            end
            
            I = flipdim(I,1);
            I = flipdim(I,2);
            if ndims(I)==4
                I = permute(I,[2 1 3 4]);
            elseif ndims(I)==3
                I = permute(I,[2 1 3]);
            end
            VDims = VDims([2 1 3]);
            
            hdr.sizeof_hdr = 348;
            hdr.db_name         = zeros(1, 18, 'uint8');
            hdr.extents         = zeros(1, 1, 'int32');
            hdr.session_error   = zeros(1, 1, 'int16');
            hdr.regular         = 'r';
            hdr.dim_info        = zeros(1, 1, 'uint8');
            hdr.dim = ones(1,8,'int16');
            hdr.dim(1) = ndims(I);
            hdr.dim(2) = size(I,1);
            hdr.dim(3) = size(I,2);
            hdr.dim(4) = size(I,3);
            hdr.dim(5) = size(I,4);
            hdr.intent_p1       = zeros(1, 1, 'single');
            hdr.intent_p2       = zeros(1, 1, 'single');
            hdr.intent_p3       = zeros(1, 1, 'single');
            hdr.intent_code     = zeros(1, 1, 'int16');
            hdr.data_type        = zeros(1, 10, 'int16');
            hdr.slice_start     = zeros(1, 1, 'int16');
            hdr.pixdim          = zeros(1, 8, 'single');
            hdr.pixdim(2) = VDims(1);
            hdr.pixdim(3) = VDims(2);
            hdr.pixdim(4) = VDims(3);
            hdr.pixdim(5:end)=1;
            hdr.vox_offset = 352;
            hdr.scl_slope       = zeros(1, 1, 'single');
            hdr.scl_inter       = zeros(1, 1, 'single');
            hdr.slice_end		= zeros(1, 1, 'int16');
            hdr.slice_code 		= zeros(1, 1, 'uint8');
            hdr.xyzt_units 		= bitor(2,16);
            hdr.cal_max 		= zeros(1, 1, 'single');
            hdr.cal_min 		= zeros(1, 1, 'single');
            hdr.slice_duration 	= zeros(1, 1, 'single');
            hdr.toffset 		= zeros(1, 1, 'single');
            hdr.glmax 			= zeros(1, 1, 'int32');
            hdr.glmax           = max(I(:));
            hdr.glmin			= zeros(1, 1, 'int32');
            hdr.descrip         = 'Created with ExploreDTI';
            hdr.descrip = [hdr.descrip  ...
                char(zeros(1, 80-length(hdr.descrip)))];
            hdr.aux_file = 'none';
            pad = zeros(1, 24-length(hdr.aux_file));
            hdr.aux_file = [hdr.aux_file  char(pad)];
            hdr.qform_code 		= zeros(1, 1, 'int16');
            hdr.sform_code  	= ones(1, 1, 'int16');
            hdr.quatern_b   	= zeros(1, 1, 'single');
            hdr.quatern_c   	= zeros(1, 1, 'single');
            hdr.quatern_d  		= zeros(1, 1, 'single');
            hdr.qoffset_x 		= zeros(1, 1, 'single');
            hdr.qoffset_y  		= zeros(1, 1, 'single');
            hdr.qoffset_z  		= zeros(1, 1, 'single');
            hdr.srow_x  		= zeros(1, 4, 'single');
            hdr.srow_y 			= zeros(1, 4, 'single');
            hdr.srow_z 			= zeros(1, 4, 'single');

            Temp = -(VDims.*[size(I,1) size(I,2) size(I,3)])/2;

            hdr.srow_x(1) = VDims(1);
            hdr.srow_x(4) = Temp(1);
            hdr.srow_y(2) = VDims(2);
            hdr.srow_y(4) = Temp(2);
            hdr.srow_z(3) = VDims(3);
            hdr.srow_z(4) = Temp(3);

            hdr.intent_name		= zeros(1, 16, 'uint8');
            hdr.magic = 'n+1     ';
            hdr.magic(4:end) = 0;

            dat = whos('I');
            switch (dat.class)
                case 'logical'
                    I = single(I);
                    hdr.datatype = 16;
                    hdr.bitpix = 32;
                case 'uint8'
                    hdr.datatype = 2;
                    hdr.bitpix = 8;
                case 'int8'
                    hdr.datatype = 256;
                    hdr.bitpix = 8;
                case 'int16'
                    hdr.datatype = 4;
                    hdr.bitpix = 16;
                case 'uint16'
                    hdr.datatype = 512;
                    hdr.bitpix = 16;
                case 'single'
                    hdr.datatype = 16;
                    hdr.bitpix = 32;
                case 'double'
                    hdr.datatype = 64;
                    hdr.bitpix = 64;
                otherwise
                    hdr.datatype = 16;
                    hdr.bitpix = 32;
            end
            
            f = fopen(file_name, 'w+');
            if f == -1
                disp('Could not write file:')
                disp(file_name)
                return;
            end
            
            fwrite(f, int32(hdr.sizeof_hdr(1)),    'int32');
            fwrite(f, uint8(hdr.data_type(1:10)),  'uint8');
            fwrite(f, uint8(hdr.db_name(1:18)),    'uint8');
            fwrite(f, int32(hdr.extents(1)),       'int32');
            fwrite(f, int16(hdr.session_error(1)), 'int16');
            fwrite(f, hdr.regular(1),              'uchar');
            fwrite(f, uint8(hdr.dim_info(1)),      'uint8');
            fwrite(f, int16(hdr.dim(1:8)),         'int16');
            fwrite(f,single(hdr.intent_p1(1)),     'float32');
            fwrite(f,single(hdr.intent_p2(1)),     'float32');
            fwrite(f,single(hdr.intent_p3(1)),     'float32');
            fwrite(f, int16(hdr.intent_code(1)),   'int16');
            fwrite(f, int16(hdr.datatype(1)),      'int16');
            fwrite(f, int16(hdr.bitpix(1)),        'int16');
            fwrite(f, int16(hdr.slice_start(1)),   'int16');
            fwrite(f,single(hdr.pixdim(1:8)),      'float32');
            fwrite(f,single(hdr.vox_offset(1)),    'float32');
            fwrite(f,single(hdr.scl_slope(1)),     'float32');
            fwrite(f,single(hdr.scl_inter(1)),     'float32');
            fwrite(f, int16(hdr.slice_end(1)),     'int16');
            fwrite(f, uint8(hdr.slice_code(1)),    'uint8');
            fwrite(f, hdr.xyzt_units,              'char');
            fwrite(f,single(hdr.cal_max(1)),       'float32');
            fwrite(f,single(hdr.cal_min(1)),       'float32');
            fwrite(f,single(hdr.slice_duration(1)),'float32');
            fwrite(f,single(hdr.toffset(1)),       'float32');
            fwrite(f, int32(hdr.glmax(1)),         'int32');
            fwrite(f, int32(hdr.glmin(1)),         'int32');
            fwrite(f, uint8(hdr.descrip(1:80)),    'uint8');
            fwrite(f, hdr.aux_file(1:24),          'uchar');
            fwrite(f, int16(hdr.qform_code(1)),    'int16');
            fwrite(f, int16(hdr.sform_code(1)),    'int16');
            fwrite(f,single(hdr.quatern_b(1)),     'float32');
            fwrite(f,single(hdr.quatern_c(1)),     'float32');
            fwrite(f,single(hdr.quatern_d(1)),     'float32');
            fwrite(f,single(hdr.qoffset_x(1)),     'float32');
            fwrite(f,single(hdr.qoffset_y(1)),     'float32');
            fwrite(f,single(hdr.qoffset_z(1)),     'float32');
            fwrite(f,single(hdr.srow_x(1:4)),      'float32');
            fwrite(f,single(hdr.srow_y(1:4)),      'float32');
            fwrite(f,single(hdr.srow_z(1:4)),      'float32');
            fwrite(f, uint8(hdr.intent_name(1:16)),'uint8');
            fwrite(f, uint8(hdr.magic(1:8)),       'uint8');
            
            switch (hdr.datatype)
                case 512
                    fwrite(f, uint16(I), 'uint16');
                case 256
                    fwrite(f, int8(I), 'int8');
                case 64
                    fwrite(f, double(I), 'float64');
                case 16
                    fwrite(f, single(I), 'float32');
                case 4
                    fwrite(f, int16(I), 'int16');
                case 2
                    fwrite(f, uint8(I), 'uint8');
            end
            
            fclose(f);
            
            clear I;
            
            if strcmpi(EXT,'.gz')
                gzip(file_name)
                delete(file_name)
            end
        end
        
        % From ExploreDTI: compute b-matrix
        function E_DTI_convert_nii_dic_2_txt_exe(bvalf,outname)
            
            bvecf = [bvalf(1:end-4) 'bvec'];
            if(nargin < 2 || isempty(outname))
                bmf = [bvalf(1:end-4) 'txt'];
            else
                bmf = outname;
            end
            
            if exist(bvecf,'file')==2
                
                try
                    gr = textread(bvecf,'%s','delimiter','\n','bufsize',10^6);
                    if length(gr)==3
                        g = [str2num(gr{1})' str2num(gr{2})' str2num(gr{3})'];
                    elseif length(str2num(gr{1}))==3
                        g = zeros(length(gr),3);
                        for i=1:length(gr)
                            g(i,:) = str2num(gr{i});
                        end
                    else
                        disp('Unknown format for gradient directions...')
                        disp(['Info: rows = ' num2str(length(gr)) '; columns = ' num2str(length(str2num(gr{1})))])
                        return;
                    end
                    
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                try
                    bval = textread(bvalf,'%f','delimiter','\n','bufsize',10^6);
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                if size(g,1)~=size(bval,1)
                    disp('Size of ''bvals'' does not match size of ''bvecs'' for files:')
                    disp([bvalf ': ' num2str(size(bval,1))])
                    disp([bvecf ': ' num2str(size(g,1))])
                    disp(' ')
                    return;
                end
                
                try
                    b = repmat(bval,[1 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) ...
                        g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                fid = fopen(bmf,'wt');
                if fid==-1
                    disp('Could not write file - check permissions...')
                else
                    for i=1:size(b,1)
                        fprintf(fid, '%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n',b(i,:));
                    end
                end
                fclose(fid);
            else
                disp('Could not find file:')
                disp(bvecf)
                disp(' ')
            end
        end       
        
        % From ExploreDTI: compute b-matrix - modified to accept different .bval /
        % .bvec names
        function E_DTI_convert_nii_dic_2_txt_exe2(bvalf,bvecf,outname)
            
            if(nargin < 2 || isempty(outname))
                bmf = [bvalf(1:end-4) 'txt'];
            else
                bmf = outname;
            end
            
            if exist(bvecf,'file')==2
                
                try
                    gr = textread(bvecf,'%s','delimiter','\n','bufsize',10^6);
                    if length(gr)==3
                        g = [str2num(gr{1})' str2num(gr{2})' str2num(gr{3})'];
                    elseif length(str2num(gr{1}))==3
                        g = zeros(length(gr),3);
                        for i=1:length(gr)
                            g(i,:) = str2num(gr{i});
                        end
                    else
                        disp('Unknown format for gradient directions...')
                        disp(['Info: rows = ' num2str(length(gr)) '; columns = ' num2str(length(str2num(gr{1})))])
                        return;
                    end
                    
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                try
                    bval = textread(bvalf,'%f','delimiter','\n','bufsize',10^6);
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                if size(g,1)~=size(bval,1)
                    disp('Size of ''bvals'' does not match size of ''bvecs'' for files:')
                    disp([bvalf ': ' num2str(size(bval,1))])
                    disp([bvecf ': ' num2str(size(g,1))])
                    disp(' ')
                    return;
                end
                
                try
                    b = repmat(bval,[1 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) ...
                        g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
                catch me
                    disp(me.message)
                    disp(' ')
                    return;
                end
                
                fid = fopen(bmf,'wt');
                if fid==-1
                    disp('Could not write file - check permissions...')
                else
                    for i=1:size(b,1)
                        fprintf(fid, '%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\t%12.8f\n',b(i,:));
                    end
                end
                fclose(fid);
            else
                disp('Could not find file:')
                disp(bvecf)
                disp(' ')
            end
        end
        
        % From ExploreDTI: masking based on image-processing
        function mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(A,tune,a)
            
            % A = 3D data set (e.g., the non-DWI).
            % tune = manual tuning factor (typical value range [0.5 1])
            % a = morphological filter size (e.g., 7)
            
            if tune==0
                mask = true(size(A));
                return;
            elseif isinf(tune)
                mask = false(size(A));
                return;
            end
            
            A = double(A);
            try
                A = imfill(A);
            catch me
                disp(me.message)
                disp('Warning: no masking is performed as it seems you')
                disp('do not have Matlab''s Image Processing toolbox installed.')
                mask = true(size(A));
                return;
            end
            
            B = A;
            A = A(:);
            A(A==0)=[];
            y = prctile(A,99);
            A(A>y)=[];
            T = median(A(:));
            
            crit = 1;
            N = 0;
            
            while crit && N<100
                N = N+1;
                T_new = (median(A(A<=T)) + median(A(A>T)))/2;
                crit = T_new ~= T;
                T = T_new;
            end
            mask = B>T*tune;
            
            [L,NUM] = bwlabeln(mask,6);
            
            if NUM==0
                mask(:)=1;
                return;
            end
            
            SL = L(L~=0);
            sd = hist(SL,1:NUM);
            [M,I] = max(sd);
            mask(L~=I)=0;
            mask = imfill(mask,'holes');
            
            se = zeros(a,a,a);
            se((a+1)/2,(a+1)/2,(a+1)/2)=1;
            se = smooth3(se,'gaussian',[a a a],1);
            se = se>=se((a+1)/2,(a+1)/2,end);
            
            mask = imerode(mask,se)>0;
            
            [L,NUM] = bwlabeln(mask,6);
            
            if NUM==0
                mask(:)=1;
                return;
            end
            
            SL = L(L~=0);
            sd = hist(SL,1:NUM);
            [M,I] = max(sd);
            mask(L~=I)=0;
            mask = imdilate(mask,se)>0;
            
            for i=1:size(mask,3)
                mask(:,:,i) = imfill(mask(:,:,i),'holes');
            end
        end
        
        % From ExploreDTI: masking based on image-processing
        function mask = E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,tune_1,tune_2,a)
            
            if tune_1==0 && tune_2==0
                mask = true(size(DWI{1}));
                return;
            end
            
            mask1 = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(DWI{1},tune_1,a);
            m2 = EDTI_Library.E_DTI_mean_DWI(DWI, NrB0);
            mask2 = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(m2,tune_2,a);
            
            mask = or(mask1,mask2);
            
            if sum(double(mask(:)))==0
                disp('Warning: mask was too small... omitting mask instead.')
                mask = true(size(DWI{1}));
            end
            
        end
        
        % From ExploreDTI: signal drift correction
        function E_DTI_signal_drift_correction(par)
            
            % par = E_DTI_settings_sdc_CL;
            
            if par.method==1
                tm = 'Linear';
            elseif par.method==2
                tm = 'Quadratic';
            elseif par.method==3
                tm = 'Cubic';
            end
            
            ttext = [tm ' fit for b-val = ' num2str(par.bvalC) ' s/mm^2'];
            
            if par.masking.do_it==1
                textm = ' (in mask)';
            else
                textm = '';
            end
            
            [DWI,VDims,suc] = EDTI_Library.E_DTI_read_nifti_file(par.f_in_nii);
            
            if suc==0
                return;
            end
            
            if ndims(DWI)~=4
                return;
            end
            
            S = size(DWI,4);
            
            try
                b = textread(par.f_in_txt);
            catch me
                disp(me.message)
                disp(par.f_in_txt)
                return;
            end
            
            if S~=size(b,1)
                disp(['Number of volumes in file ''' par.f_in_nii ''': ' num2str(S)])
                disp(['Number of rows in b-matrix file ''' par.f_in_txt ''': ' num2str(size(b,1))])
                disp('These values should be equal.')
                return;
            end
            
            bval = sum(b(:,[1 4 6]),2);
            
            F = find(abs(bval-par.bvalC)<=par.bv_thresh);
            G = find(abs(bval-par.bvalC)>par.bv_thresh);
            
            if isempty(F)
                disp(['Could not find data with b-values of ' num2str(par.bvalC,4) 's/mm^2'])
                disp('Processing unsuccessful for file:')
                disp(par.f_in_nii)
                return;
            end
            
            if par.method >= length(F)
                disp('Not enough data points to perform the fit.')
                disp('Processing unsuccessful for file:')
                disp(par.f_in_nii)
                disp('Adjust settings.')
                return;
            end
            
            if isempty(G)
                disp('Warning: fit may not be optimal for file:')
                disp(par.f_in_nii)
                disp('Adjust settings.')
            end
            
            for i=1:length(F)
                b0{i} = DWI(:,:,:,F(i));
                if par.masking.do_it==1
                    maskb0{i} = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(b0{i},par.masking.p2,par.masking.p1);
                else
                    maskb0{i} = true(size(b0{i}));
                end
            end
            % mask_all_b0s = false(size(DWI(:,:,:,1)));
            % for i=1:length(F)
            %     mask_all_b0s = or(mask_all_b0s, maskb0{i});
            % end
            % for i=1:length(F)
            %     maskb0{i} = mask_all_b0s;
            % end
            for i=1:length(F)
                b0s(i,1) = mean(b0{i}(maskb0{i}));
            end
            
            
            for i=1:length(G)
                dw{i} = DWI(:,:,:,G(i));
                if par.masking.do_it==1
                    maskdwi{i} = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(dw{i},par.masking.p2,par.masking.p1);
                else
                    maskdwi{i} = true(size(dw{i}));
                end
            end
            % mask_all_dws = false(size(DWI(:,:,:,1)));
            % for i=1:length(G)
            %     mask_all_dws = or(mask_all_dws, maskdwi{i});
            % end
            % for i=1:length(G)
            %     maskdwi{i} = mask_all_dws;
            % end
            for i=1:length(G)
                dws(i,1) = mean(dw{i}(maskdwi{i}));
            end
            
            warning off all
            if par.method==1
                fitr = fit(F, b0s, 'a*x+b');
                fac = (1:S).*fitr.a + fitr.b;
                normf = fitr.b;
            elseif par.method==2
                fitr = fit(F, b0s, 'a*x^2+b*x+c');
                fac = ((1:S).^2.*fitr.a) + (1:S).*fitr.b + fitr.c;
                normf = fitr.c;
            elseif par.method==3
                fitr = fit(F, b0s, 'a*x^3+b*x^2+c*x+d');
                fac = ((1:S).^3.*fitr.a) + (1:S).^2.*fitr.b + (1:S).*fitr.c + fitr.d;
                normf = fitr.d;
            end
            warning on all
            
            
            corr_fac = (normf./fac);
            
            perc_change = ((fac-normf)/abs(normf))*100;
            [mp,Ip] = max(abs(perc_change));
            perc_change = perc_change(Ip);
            
            DWIc = single(DWI);
            
            for i=1:S
                DWIc(:,:,:,i) = DWI(:,:,:,i)*corr_fac(i);
            end
            
            for i=1:length(F)
                
                b0c{i} = DWIc(:,:,:,F(i));
                
                if par.masking.do_it==1
                    b0cs(i,1) = mean(b0c{i}(maskb0{i}));
                else
                    b0cs(i,1) = mean(b0c{i}(:));
                end
                
            end
            
            for i=1:length(G)
                
                dwc{i} = DWIc(:,:,:,G(i));
                
                if par.masking.do_it==1
                    dwcs(i,1) = mean(dwc{i}(maskdwi{i}));
                else
                    dwcs(i,1) = mean(dwc{i}(:));
                end
                
            end
            
            M = 0;
            M = max(max(b0s),M);
            M = max(max(b0cs),M);
            if exist('dws','var')==1
                M = max(max(dws),M);
                M = max(max(dwcs),M);
            end
            
            m = inf;
            m = min(min(b0s),m);
            m = min(min(b0cs),m);
            if exist('dws','var')==1
                m = min(min(dws),m);
                m = min(min(dwcs),m);
            end
            
            dif = (M-m)/10;
            m = m-dif;
            M = M+dif;
            
            if perc_change<0
                text_pc = [' (' num2str(abs(perc_change),2) '% signal decrease)'];
            else
                text_pc = [' (' num2str(abs(perc_change),2) '% signal increase)'];
            end
            
            if par.show_summ_plot>0
                
                h_k = figure('Name','Drift correction','NumberTitle' ,'off','Color','w','visible','off');
                %     EDTI_Library.E_DTI_Set_Fig_Icon(h_k);
                plot(F,b0s,'or','MarkerSize',4,'Markerfacecolor',[0.8 0.8 0.8])
                hold on;
                if exist('dws','var')==1
                    plot(G,dws,'dr','MarkerSize',6,'Markerfacecolor',[0.8 0.8 0.8])
                end
                plot(F,b0cs,'ob','MarkerSize',4,'Markerfacecolor',[0.8 0.8 0.8])
                if exist('dws','var')==1
                    plot(G,dwcs,'db','MarkerSize',6,'Markerfacecolor',[0.8 0.8 0.8])
                end
                plot(1:S,fac,'-k')
                plot([-2 S+3], [normf normf],...
                    'linestyle','--','color',[0.5 0.5 0.5]); xlim([-1 S+2]);
                ax = axis;
                axis([ax(1) ax(2) m M])
                xlabel('Diffusion-weighted volume number')
                ylabel(['Mean signal per volume' textm])
                title([ttext text_pc]);
                if exist('dws','var')==1
                    legend('Uncorrected (used for fit)', 'Uncorrected (not used for fit)',...
                        'Corrected (from fit)','Corrected (not from fit)','Fit','reference','Location','west')
                else
                    legend('Uncorrected (used for fit)',...
                        'Corrected (from fit)','Fit','reference','Location','west')
                end
                set(h_k, 'InvertHardCopy', 'off');
                print(h_k,'-dpng','-r300',[par.f_out_nii(1:end-3) 'png']);
                if par.show_summ_plot==2
                    set(h_k,'visible','on')
                    savefig(h_k,[par.f_out_nii(1:end-3) 'fig'])
                end
                close(h_k)
            end
            
%             if max(DWIc(:))<intmax('int16')
%                 clear DWI;
%                 DWIc = int16(round(DWIc));
%             end
            
            EDTI_Library.E_DTI_write_nifti_file(DWIc,VDims,par.f_out_nii);
            copyfile(par.f_in_txt,par.f_out_txt);
            
        end
        
        % From ExploreDTI: flip / permute the NIFTI image matrix. The rotation
        % matrices are discarded
        function E_DTI_flip_permute_nii_file_exe(f,p,f_out)
            
            try
                [A, VDims] = EDTI_Library.E_DTI_read_nifti_file(f);
            catch
                disp('Error (memory/permission issues?) occured while trying to load file:')
                disp(f)
                return;
            end
            
            if ndims(A)~=4 && ndims(A)~=3
                disp(['File ''' f ''' is not 3D or 4D...'])
                return;
            end
            
            A = permute(A,[p.permute 4]);
            VDims = VDims(p.permute);
            
            for i=1:3
                if p.flip(i)==1
                    A = flipdim(A,i);
                end
            end
            
            if ~isempty(p.force_voxel_size)
                VDims = p.force_voxel_size;
            end
            
            try
                EDTI_Library.E_DTI_write_nifti_file(A,VDims,f_out);
            catch
                disp('Error (permission?) occured while writing file:')
                disp(f_out)
                return;
            end
        end
        
        % From ExploreDTI: flip / permute the b-matrix as specified
        function [g, b] = E_DTI_flip_g_b(g,b,perm,flip)
            
            if perm==2
                g(:,[1 2]) = g(:,[2 1]);
                b(:,[1 2 3 4 5 6]) = b(:,[4 2 5 1 3 6]);
            elseif perm==3
                g(:,[1 3]) = g(:,[3 1]);
                b(:,[1 2 3 4 5 6]) = b(:,[6 5 3 4 2 1]);
            elseif perm==4
                g(:,[2 3]) = g(:,[3 2]);
                b(:,[1 2 3 4 5 6]) = b(:,[1 3 2 6 5 4]);
            elseif perm==5
                g(:,[1 2 3]) = g(:,[2 3 1]);
                b(:,[1 2 3 4 5 6]) = b(:,[4 3 5 6 2 1]);
            elseif perm==6
                g(:,[1 2 3]) = g(:,[3 1 2]);
                b(:,[1 2 3 4 5 6]) = b(:,[6 3 5 1 2 4]);
            end
            
            if flip==2
                g(:,1) = -g(:,1);
                b(:,[2 3]) = -b(:,[2 3]);
            elseif flip==3
                g(:,2) = -g(:,2);
                b(:,[2 5]) = -b(:,[2 5]);
            elseif flip==4
                g(:,3) = -g(:,3);
                b(:,[3 5]) = -b(:,[3 5]);
            end
        end
        
        % From ExploreDTI: get the gradients from the b-matrix
        function [bval, grad] =  E_DTI_GetGradientsandBval_SC(bmat, NrB0)
            
            bmat = double(bmat);
            bval = (sum(bmat(NrB0+1:end,[1 4 6]),2));
            
            % x~=0
            BSign_x = sign(sign(bmat(:,1:3)) + 0.0001);   % Adding 0.0001 avoids getting zeros here
            
            % y~=0
            BSign_y = sign(sign(bmat(:,[2 4 5])) + 0.0001);   % Adding 0.0001 avoids getting zeros here
            
            Bo = bmat(:,1)==0;
            
            BSign = BSign_x;
            BSign([Bo Bo Bo]) = BSign_y([Bo Bo Bo]);
            
            grad = BSign .* sqrt(bmat(:,[1 4 6]));
            grad_n = sqrt(sum(grad.*grad,2));
            grad = grad./[grad_n grad_n grad_n];
            
            if ~isreal(grad)
                %     disp('Warning: incorrect values of the B-matrix were encountered!')
                grad = real(grad);
                grad_n = sqrt(sum(grad.*grad,2));
                grad = grad./[grad_n grad_n grad_n];
            end
            
            grad(isnan(grad))=0;
            
            grad = grad(NrB0+1:end,:);
        end
        
        % From ExploreDTI: DTI scalars and first eigenvector
        function [FEFA, FA, FE, SE, eigval, Bm] = E_DTI_eigensystem_analytic(DT)
            
            mask = ~isnan(DT{1});
            
            a1 = DT{1}(mask);
            a2 = DT{2}(mask);
            a3 = DT{3}(mask);
            a4 = DT{4}(mask);
            a5 = DT{5}(mask);
            a6 = DT{6}(mask);
            
            FE = nan([size(mask) 3]);
            SE = FE;
            eigval = FE;
            Bm = false([size(mask) 7]);
            
            I1 = a1 + a4 + a6;
            I2 = a1.*a4 + a1.*a6 + a4.*a6;
            I2 = I2 - (a2.*a2 + a3.*a3 + a5.*a5);
            I3 = a1.*a4.*a6 + 2*a2.*a3.*a5;
            I3 = I3 - (a6.*a2.*a2 + a4.*a3.*a3 + a1.*a5.*a5);
            
            v = (I1.*I1)/9 - I2/3;
            s = (I1.*I1.*I1)/27 - (I1.*I2)/6 + I3/2;
            
            vt = v;
            vt(vt<=0)=1;
            
            temp = ((s./vt).*(1./sqrt(vt)));
            
            temp(temp<-1)=-1;
            temp(temp>1)=1;
            
            phi = acos(temp)/3;
            
            E{1} = I1/3 + 2*sqrt(v).*cos(phi);
            E{2} = I1/3 - 2*sqrt(v).*cos(pi/3 + phi);
            E{3} = I1/3 - 2*sqrt(v).*cos(pi/3 - phi);
            
            
            for i=1:3
                
                A = a1-E{i};
                B = a4-E{i};
                C = a6-E{i};
                
                Vx{i} = (a2.*a5 - B.*a3).*(a3.*a5 - C.*a2);
                Vy{i} = (a3.*a5 - C.*a2).*(a3.*a2 - A.*a5);
                Vz{i} = (a2.*a5 - B.*a3).*(a3.*a2 - A.*a5);
                
            end
            
            
            Q{1} = I3>0;
            Q{2} = a1>=0;
            Q{3} = a4>=0;
            Q{4} = a6>=0;
            Q{5} = (a1.*a4 - a2.*a2)>=0;
            Q{6} = (a1.*a6 - a3.*a3)>=0;
            Q{7} = (a4.*a6 - a5.*a5)>=0;
            
            dummy = nan(size(mask));
            
            for i=1:3
                dummy(mask) = E{i};
                eigval(:,:,:,i) = dummy;
            end
            
            NV1 = sqrt(Vx{1}.^2 + Vy{1}.^2 + Vz{1}.^2);
            MaskV = NV1==0;
            NV1(MaskV) = 1;
            Vx{1}(MaskV) = 1;
            Vy{1}(MaskV) = 0;
            Vz{1}(MaskV) = 0;
            
            NV2 = sqrt(Vx{2}.^2 + Vy{2}.^2 + Vz{2}.^2);
            MaskV = NV2==0;
            NV2(MaskV) = 1;
            Vx{2}(MaskV) = 1;
            Vy{2}(MaskV) = 0;
            Vz{2}(MaskV) = 0;
            
            Vx{1} = Vx{1}./NV1;
            Vy{1} = Vy{1}./NV1;
            Vz{1} = Vz{1}./NV1;
            Vx{2} = Vx{2}./NV2;
            Vy{2} = Vy{2}./NV2;
            Vz{2} = Vz{2}./NV2;
            
            
            dummy(mask) = Vx{1};
            FE(:,:,:,1) = dummy;
            dummy(mask) = Vy{1};
            FE(:,:,:,2) = dummy;
            dummy(mask) = Vz{1};
            FE(:,:,:,3) = dummy;
            
            dummy(mask) = Vx{2};
            SE(:,:,:,1) = dummy;
            dummy(mask) = Vy{2};
            SE(:,:,:,2) = dummy;
            dummy(mask) = Vz{2};
            SE(:,:,:,3) = dummy;
            
            
            FA = EDTI_Library.FrAn_calc(eigval(:,:,:,1),eigval(:,:,:,2),eigval(:,:,:,3));
            
            FEFA = nan([size(FA) 3]);
            
            for i=1:3
                dummy = FE(:,:,:,i);
                dummy(mask) = FA(mask).*abs(dummy(mask));
                FEFA(:,:,:,i) = dummy;
            end
            
            FEFA(FEFA<0)=0;
            FEFA(FEFA>1)=1;
            FA(mask) = FA(mask)*sqrt(3);
            FA(FA>sqrt(3))=sqrt(3);
            FA(FA<0)=0;
            
            dummy = false(size(mask));
            
            for i=1:length(Q)
                dummy(mask) = ~Q{i};
                Bm(:,:,:,i) = dummy;
            end
        end
        
        % From ExploreDTI: helper function (cosine angle)
        function a = angle2(v1, v2)
            % v1 = normalize(v1);
            % v2 = normalize(v2);
            a = 180/pi*real(acos(abs(sum(v1.*v2,1))));
        end
        
        % From ExploreDTI: Fractional Anisotropy
        function fa = FrAn_calc(L1,L2,L3)
            
            mask = ~isnan(L1);
            L1 = L1(mask);
            L2 = L2(mask);
            L3 = L3(mask);
            
            D1 = L1-L2;
            D1 = D1.^2;
            
            D2 = L2-L3;
            D2 = D2.^2;
            
            D1 = D1 + D2;
            
            D2 = L3-L1;
            D2 = D2.^2;
            
            D1 = D1 + D2;
            D1 = sqrt(D1);
            
            clear D2;
            
            L1 = L1.^2;
            L2 = L2.^2;
            L3 = L3.^2;
            
            L1 = L1 + L2;
            L1 = L1 + L3;
            
            clear L2 L3;
            
            L1 = 2*L1;
            L1 = sqrt(L1);
            fa = repmat(double(nan),size(mask));
            
            warning off all
            dummy = D1./L1;
            dummy(isnan(dummy))=0;
            dummy(isinf(dummy))=0;
            fa(mask) = dummy;
            warning on all
            
        end
        
        % From ExploreDTI: Compute outliers using the residuals
        function [outlier, chi_sq, chi_sq_iqr, residuals] = E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par)
            
            outlier = false([size(DWIB0) size(b,1)]);
            chi_sq = [];
            chi_sq_iqr = [];
            
            mask = ~isnan(DT{1});
            MDims = size(mask);
            
            Vm = sum(mask(:));
            
            X = repmat(double(0),[6 Vm]);
            
            for k=1:6
                X(k,:) = DT{k}(mask(:));
            end
            clear DT;
            
            DWI_m = repmat(double(0),[size(b,1) Vm]);
            
            for k=1:length(DWI)
                DWI_m(k,:) = double(DWI{k}(mask(:)));
            end
            clear DWI;
            
            B0_est = DWIB0(mask(:))';
%             clear DWIB0;
            
            res = -double(b)*X;
            
            clear X;
            
            res = exp(res);
            res(isinf(res)) = max(res(~isinf(res)));
            
            for i=1:size(res,1)
                res(i,:) = res(i,:).*B0_est;
            end
            
            res = res - DWI_m;
            residuals = zeros(size(res,1),size(DWIB0,1)*size(DWIB0,2)*size(DWIB0,3));
            residuals(:,mask) = res;
            residuals = reshape(residuals',size(DWIB0,1),size(DWIB0,2),size(DWIB0,3),size(res,1));
            res = abs(res);
            clear DWI_m;
            
            m = prctile(res,[25 50 75],1);
            chi_sq = repmat(double(nan),MDims);
            chi_sq(mask) = m(2,:);
            
            chi_sq_iqr = repmat(double(nan),MDims);
            chi_sq_iqr(mask) = m(3,:)-m(1,:);
            
            for i=1:size(outlier,4)
                
                dummy = false(MDims);
                dummy(mask) = res(i,:)' > m(3,:)' + par.RE.kappa*chi_sq_iqr(mask);
                outlier(:,:,:,i)=dummy;
                
            end
        end
        
        % From ExploreDTI: Compute outliers using the residuals (DKI version)
        function [outlier, chi_sq, chi_sq_iqr,residuals] = E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT)
            
            b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
            
            outlier = false([size(DWIB0) size(b,1)]);
            chi_sq = [];
            chi_sq_iqr = [];
            
            mask = ~isnan(DT{1});
            MDims = size(mask);
            
            Vm = sum(mask(:));
            
            X = repmat(single(0),[6 Vm]);
            
            for k=1:6
%                 X(k,:) = DT{k}(mask(:));
            end
            clear DT;
            
            Y = repmat(single(0),[15 Vm]);
            
            for k=1:15
                Y(k,:) = KT{k}(mask(:));
            end
            clear KT;
            
            DWI_m = repmat(single(0),[size(b,1) Vm]);
            
            for k=1:length(DWI)
                DWI_m(k,:) = single(DWI{k}(mask(:)));
            end
            clear DWI;
            
            B0_est = DWIB0(mask(:))';
%             clear DWIB0;
            
            res = -single(b)*X + (b_kurt*Y).*repmat((X(1,:)+X(4,:)+X(6,:)).^2,[size(b_kurt,1) 1]);
            clear X Y;
            
            res = exp(res);
            
            for i=1:size(res,1)
                res(i,:) = res(i,:).*B0_est;
            end
            
            res = res - DWI_m;
            residuals = zeros(size(res,1),size(DWIB0,1)*size(DWIB0,2)*size(DWIB0,3));
            residuals(:,mask) = res;
            residuals = reshape(residuals',size(DWIB0,1),size(DWIB0,2),size(DWIB0,3),size(res,1));
            res = abs(res);
            clear DWI_m;
            
            m = prctile(res,[25 50 75],1);
            
            chi_sq = repmat(single(nan),MDims);
            chi_sq(mask) = m(2,:);
            
            chi_sq_iqr = repmat(single(nan),MDims);
            chi_sq_iqr(mask) = m(3,:)-m(1,:);
            
            for i=1:size(outlier,4)
                
                dummy = false(MDims);
                dummy(mask) = res(i,:)' > m(3,:)' + par.RE.kappa*chi_sq_iqr(mask);
                outlier(:,:,:,i)=dummy;
                
            end
        end
        
        % From ExploreDTI: Weighted least squares
        function X = E_DTI_WLLS_WW(S0,B)
            
            warning off all
            
            
            crit = 1;
            iter_max = 10;
            cn=0;
            con = 10^-4;
            
            Log_S0 = log(S0);
            
            W = ones(length(S0),1);
            
            
            
            X = ones(size(B,2),1);
            
            
            while crit==1
                
                cn = cn+1;
                W_p = W;
                X_p = X;
                w = diag(W_p.^2);
                X = (B'*w*B)\(B'*w)*Log_S0;
                W = exp(B*X);
                
                if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn==iter_max
                    crit=0;
                end
                
            end
            
            if any(isnan(X)) || any(isinf(X))
                
                X = B\S0;
                
            end
            
            
            % disp(num2str(cn))
            
            warning on all
        end
        
        % From ExploreDTI: REKINDLE
        function [X, is_outlier] = E_DTI_robust_Linear_fit_new(S0,B,kappa,iter_max,con)
            
            warning off all
            
            % Initialization
            crit_all = 1;
            cn_all = 0;
            Max_inter_in = 5;
            
            % Step 1: Initial LLS fit
            w = ones(size(S0));
            LogS0 = log(S0);
            W = diag(w);
            X = (B'*W*B)\(B'*W)*LogS0;
            
            while crit_all==1
                
                % Initialization
                cn_all = cn_all + 1;
                X_p_all = X;
                cn = 0;
                crit = 1;
                
                % Step 2: Compute a robust estimate for homoscedastic regression using IRLS.
                while crit==1
                    
                    % Initialization
                    cn = cn+1;
                    X_p = X;
                    
                    fit = B*X;
                    % a. Calculate the residuals e in the linear domain
                    res = (LogS0-fit);
                    if all(res==0)
                        break;
                    end
                    % b. Obtain an estimate of the dispersion of the residuals by
                    % calculating the median absolute deviation (MAD).
                    C = 1.4826*median(abs(res-median(res))) ;
                    % c. Recompute the weights according to Eq. [13].
                    w = 1./(1 + (res/C).^2).^2;
                    W = diag(w);
                    % d. Perform WLLS fit with new weights
                    X = (B'*W*B)\(B'*W)*LogS0;
                    % e. Check convergence
                    if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == Max_inter_in
                        crit=0;
                    end
                    
                end
                
                % Step 3: Transform variables for heteroscedasticity
                fit = B*X;
                LogS0_2 = LogS0./(exp(-fit));
                B2 = B./repmat((exp(-fit)),[1,size(B,2)]);
                
                % Initialization
                crit = 1;
                cn = 0;
                
                % Step 4: Initial LLS fit in * domain
                w = ones(size(S0));
                W = diag(w);
                X = (B'*W*B)\(B'*W)*LogS0;
                
                % Step 5: Compute a robust estimate for homoscedastic regression using IRLS.
                while crit==1
                    
                    % Initialization
                    cn = cn+1;
                    X_p = X;
                    
                    fit = B2*X;
                    % a. Calculate the residuals e* in the linear domain
                    res = (LogS0_2-fit);
                    if all(res==0)
                        break;
                    end
                    % b. Obtain an estimate of the dispersion of the residuals by
                    % calculating the median absolute deviation (MAD).
                    C = 1.4826*median(abs(res-median(res))) ;
                    % c. Recompute the weights according to Eq. [13].
                    w = 1./(1 + (res/C).^2).^2;
                    W = diag(w);
                    
                    % d. Perform WLLS fit with new weights
                    X = (B2'*W*B2)\(B2'*W)*LogS0_2;
                    % e. Check convergence
                    if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == Max_inter_in
                        crit=0;
                    end
                    
                end
                
                %  Step 6: Check convergence overall loop
                if all(abs(X-X_p_all) <= con*max(abs(X),abs(X_p_all))) || cn_all == iter_max
                    crit_all=0;
                end
            end
            
            % Step 7: Identify and exclude outliers
            fit = B2*X; % In the first iteration, this is the first fit
            res = (LogS0_2-fit);
            
            C = 1.4826*median(abs(res-median(res))) ;
            IND = (abs(res)<= kappa*C);
            
            % Step 8: Final fit
            X = EDTI_Library.E_DTI_WLLS_WW(S0(IND),B(IND,:));
            
            is_outlier = ~IND;
            
            warning on all
        end
        
        % From ExploreDTI: Weighted NLLS
        function X = E_DTI_quick_par(options,dwi,b_final)
            
            warning off all
            covar = diag(dwi.^2);
            
            X = (b_final'*covar*b_final)\(b_final'*covar)*log(dwi);
            
            try
                X = nlinfit(b_final,dwi,@(p,x)exp(x*p),X,options);
            catch
            end
            warning on all
        end
        
        % From ExploreDTI: Min-Max of B0s
        function [M, M1, M2] = E_DTI_phpv(DWI, NrB0)
            
            if ~iscell(DWI)
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            end
            
            pp = size(DWI{1});
            M1 = repmat(single(0),pp);
            M2 = repmat(single(0),pp);
            
            dum = repmat(single(0),[pp(1) pp(2) NrB0]);
            
            for j=1:pp(3)
                
                for i=1:NrB0
                    
                    dum(:,:,i) = DWI{i}(:,:,j);
                    
                end
                
                M1(:,:,j) = min(dum,[],3);
                
            end
            
            dum = repmat(single(0),[pp(1) pp(2) length(DWI)-NrB0]);
            
            for j=1:pp(3)
                
                for i=NrB0+1:length(DWI)
                    
                    dum(:,:,i) = DWI{i}(:,:,j);
                    
                end
                
                M2(:,:,j) = max(dum,[],3);
                
            end
            
            M = single(M1<=M2);
        end
        
        % From ExploreDTI: Max of non B0s
        function M_B0 = E_DTI_max_B0s(DWI)
            
            if ~iscell(DWI)
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            end
            
            pp = size(DWI{1});
            M_B0 = repmat(single(0),pp);
            
            
            dum = repmat(single(0),[pp(1) pp(2) length(DWI)]);
            
            for j=1:pp(3)
                
                for i=1:length(DWI)
                    
                    dum(:,:,i) = DWI{i}(:,:,j);
                    
                end
                
                M_B0(:,:,j) = max(dum,[],3);
                
            end
            
            M_B0 = single(M_B0)+1;
        end
        
        % From ExploreDTI: Remove Physically Implausible Signals (PIS)
        function B0s = E_DTI_Clean_up_B0s_2(DWI, mask, NrB0)
            
            [M, M1, M2] = EDTI_Library.E_DTI_phpv(DWI, NrB0);
            
            B0s = DWI(1:NrB0);
            M_B0 = EDTI_Library.E_DTI_max_B0s(DWI(NrB0+1:end));
            clear DWI;
            
            nan_Mask = M;
            nan_Mask(~mask)=0;
            
            I=size(M,1);
            J=size(M,2);
            K=size(M,3);
            
            for i=1:I
                for j=1:J
                    for k=1:K
                        if nan_Mask(i,j,k)
                            for u=1:length(B0s)
                                if B0s{u}(i,j,k) <= M_B0(i,j,k)
                                    B0s{u}(i,j,k) = M_B0(i,j,k);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % From ExploreDTI: The actual tensor fit
        function [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par,NrB0)
            
            if par.clean_up_PIS==1
                if NrB0>0
                    B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, mask, NrB0);
                    DWI(1:NrB0)=B0s;
                    clear B0s;
                end
            end
            
            if par.TE==1
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                b2 = [ones(size(b,1),1) -b];
                X = b2\log(DWI_m);
                
                DWIB0 = nan(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
                
            elseif par.TE==20
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                b2 = [ones(size(b,1),1) -b];
                
                X = zeros(7,Vm);
                
                parfor i=1:Vm
                    
                    dwi = DWI_m(:,i);
                    
                    covar = diag(dwi.^2);
                    
                    X(:,i) = inv(b2'*covar*b2)*(b2'*covar)*log(dwi);
                    
                end
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                DWIB0 = nan(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
                
            elseif par.TE==2
                
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                b2 = [ones(size(b,1),1) -b];
                
                X = zeros(7,Vm);
                
                parfor i=1:Vm
                    
                    dwi = DWI_m(:,i);
                    X(:,i) = EDTI_Library.E_DTI_WLLS_WW(dwi,b2);
                    
                end
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                DWIB0 = nan(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
                
            elseif par.TE==3 || (par.TE==4 && par.ROBUST_option==2)
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                b2 = [ones(size(b,1),1) -b];
                
                X = zeros(7,Vm);
                
                options = statset('nlinfit');
                if par.TE==4
                    options.Robust = 'on';
                end
                
                parfor i=1:Vm
                    
                    X(:,i) = EDTI_Library.E_DTI_quick_par(options,DWI_m(:,i),b2);
                    
                end
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                DWIB0 = nan(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
                
            elseif par.TE==4 && par.ROBUST_option==1
                
                con = par.RE.rel_convergence;
                iter_max = par.RE.max_iter;
                kappa = par.RE.kappa;
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                b2 = [ones(size(b,1),1) -b];
                
                X = zeros(7,Vm);
                outlier_m = false([length(DWI) Vm]);
                outlier = false([size(mask) length(DWI)]);
                
                parfor i=1:Vm
                    %         X(:,i) = EDTI_Library.E_DTI_robust_Linear_fit_old(DWI_m(:,i),b2,kappa,iter_max,con);
                    [X(:,i), outlier_m(:,i)] = EDTI_Library.E_DTI_robust_Linear_fit_new(DWI_m(:,i),b2,kappa,iter_max,con);
                end
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                dummy = false(size(DWI{1}));
                
                for i=1:length(DWI)
                    dummy(mask) = outlier_m(i,:);
                    outlier(:,:,:,i) = dummy;
                end
                
                DWIB0 = nan(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [dummy_stuff, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
                
            end
            
        end
        
        % From ExploreDTI: kurtosis constraints
        function mask = E_DTI_constrain_KT_mask(KT,constr1,constr2)
            
            % constr1 = [0 3];
            % constr2  = [-3 3];
            
            mask(:,:,:,1) = or(KT{1}<constr1(1),KT{1}>constr1(2));
            mask(:,:,:,2) = or(KT{2}<constr1(1),KT{2}>constr1(2));
            mask(:,:,:,3) = or(KT{3}<constr1(1),KT{3}>constr1(2));
            
            mask(:,:,:,4) = or(KT{10}<constr1(1),KT{10}>constr1(2));
            mask(:,:,:,5) = or(KT{11}<constr1(1),KT{11}>constr1(2));
            mask(:,:,:,6) = or(KT{12}<constr1(1),KT{12}>constr1(2));
            
            mask(:,:,:,7) = or(KT{4}<constr2(1),KT{4}>constr2(2));
            mask(:,:,:,8) = or(KT{5}<constr2(1),KT{5}>constr2(2));
            mask(:,:,:,9) = or(KT{6}<constr2(1),KT{6}>constr2(2));
            mask(:,:,:,10) = or(KT{7}<constr2(1),KT{7}>constr2(2));
            mask(:,:,:,11) = or(KT{8}<constr2(1),KT{8}>constr2(2));
            mask(:,:,:,12) = or(KT{9}<constr2(1),KT{9}>constr2(2));
            mask(:,:,:,13) = or(KT{13}<constr2(1),KT{13}>constr2(2));
            mask(:,:,:,14) = or(KT{14}<constr2(1),KT{14}>constr2(2));
            mask(:,:,:,15) = or(KT{15}<constr2(1),KT{15}>constr2(2));
            
            mask = any(mask,4);
        end
        
        % From ExploreDTI: The actual kurtosis fit
        function [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0)
            
            if par.clean_up_PIS==1
                if NrB0>0
                    B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, mask, NrB0);
                    DWI(1:NrB0)=B0s;
                    clear B0s;
                end
            end
            
            if par.TE==1
                
                b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                % [W1111 W2222 W3333 W1112 W1113 W2221 W2223 W3331 W3332 W1122 W1133 W2233 W1123 W2213 W3312]
                
                b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
                
                b_final = [ones(length(DWI),1) -b b_kurt];
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                X = b_final\log(DWI_m);
                
                DWIB0 = zeros(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                for i=7:21;
                    KT{i-6} = dummy;
                    KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
                
            elseif par.TE==20
                
                b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
                
                b_final = [ones(length(DWI),1) -b b_kurt];
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                X = zeros(22,Vm);
                
                parfor i=1:Vm
                    
                    dwi = DWI_m(:,i);
                    
                    covar = diag(dwi.^2);
                    
                    X(:,i) = (b_final'*covar*b_final)\(b_final'*covar)*log(dwi);
                    
                end
                
                DWIB0 = zeros(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                for i=7:21;
                    KT{i-6} = dummy;
                    KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
                
            elseif par.TE==2
                
                
                b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
                
                b_final = [ones(length(DWI),1) -b b_kurt];
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                X = zeros(22,Vm);
                
                parfor i=1:Vm
                    
                    X(:,i) = EDTI_Library.E_DTI_WLLS_WW(DWI_m(:,i),b_final);
                    
                end
                
                DWIB0 = zeros(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                for i=7:21;
                    KT{i-6} = dummy;
                    KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
                
            elseif par.TE==3 || (par.TE==4 && par.ROBUST_option==2)
                
                b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
                
                b_final = [ones(length(DWI),1) -b b_kurt];
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                X = zeros(22,Vm);
                
                options = statset('nlinfit');
                if par.TE==4
                    options.Robust = 'on';
                end
                
                parfor i=1:Vm
                    
                    X(:,i) = EDTI_Library.E_DTI_quick_par(options,DWI_m(:,i),b_final);
                    
                end
                
                DWIB0 = zeros(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                for i=7:21;
                    KT{i-6} = dummy;
                    KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
                
            elseif par.TE==4 && par.ROBUST_option==1
                
                con = par.RE.rel_convergence;
                iter_max = par.RE.max_iter;
                kappa = par.RE.kappa;
                
                b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
                
                b_final = [ones(length(DWI),1) -b b_kurt];
                
                Vm = sum(mask(:));
                DWI_m = repmat(double(0),[length(DWI) Vm]);
                
                for i=1:length(DWI)
                    DWI_m(i,:) = double(DWI{i}(mask));
                end
                
                DWI_m(DWI_m<=0)=1;
                X = zeros(22,Vm);
                
                outlier_m = false([length(DWI) Vm]);
                outlier = false([size(mask) length(DWI)]);
                
                parfor i=1:Vm
                    
                    %         X(:,i) = EDTI_Library.E_DTI_robust_Linear_fit_old(DWI_m(:,i),b_final,kappa,iter_max,con);
                    [X(:,i), outlier_m(:,i)] = EDTI_Library.E_DTI_robust_Linear_fit_new(DWI_m(:,i),b_final,kappa,iter_max,con);
                    
                end
                
                DWIB0 = zeros(size(DWI{1}));
                DWIB0(mask) = exp(X(1,:));
                
                dummy = nan(size(DWI{1}));
                
                for i=1:6
                    DT{i} = dummy;
                    DT{i}(mask) = X(i+1,:);
                end
                
                for i=7:21;
                    KT{i-6} = dummy;
                    KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
                end
                
                dummy = false(size(DWI{1}));
                
                for i=1:length(DWI)
                    dummy(mask) = outlier_m(i,:);
                    outlier(:,:,:,i) = dummy;
                end
                
                dwib0 = DWIB0(mask);
                y = max(dwib0(~isinf(dwib0)));
                dwib0(dwib0>y)=y;
                DWIB0(mask) = dwib0;
                
                [~, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
                
            end
            
        end
        
        % From ExploreDTI: Kurtosis fit wrapper
        function [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par,NrB0,VDims)
            
            [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0);
            
            if par.DKI_constraints.do_it==1
                
                mask_Fin = mask;
                
                constr1 = par.DKI_constraints.constr1;
                constr2 = par.DKI_constraints.constr2;
                
                FWHM = max(VDims);
                sigma_sm = FWHM/(2*sqrt(2*log(2)));
                
                mask = EDTI_Library.E_DTI_constrain_KT_mask(KT,constr1,constr2);
                
                crit = sum(single(mask(:)))>0;
                cn = 0;
                
                while crit
                    
                    cn=cn+1;
                    
                    for t=1:length(DWI)
                        DWI{t} = EDTI_Library.E_DTI_smooth3_ani_voxel_size(single(DWI{t}), 'gaussian', [3 3 3], sigma_sm, VDims);
                    end
                    
                    [DT_, DWIB0_, KT_, outlier_, chi_sq_, chi_sq_iqr_] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0);
                    
                    for i=1:length(KT)
                        KT{i}(mask) = KT_{i}(mask);
                    end
                    
                    for i=1:length(DT)
                        DT{i}(mask) = DT_{i}(mask);
                    end
                    
                    DWIB0(mask) = DWIB0_(mask);
                    chi_sq(mask) = chi_sq_(mask);
                    chi_sq_iqr(mask) = chi_sq_iqr_(mask);
                    outlier(repmat(mask,[1 1 1 size(outlier,4)])) = outlier_(repmat(mask,[1 1 1 size(outlier,4)]));
                    
                    mask = EDTI_Library.E_DTI_constrain_KT_mask(KT,constr1,constr2);
                    
                    
                    if cn==10 || sum(single(mask(:)))==0
                        crit=0;
                    end
                    
                end
                
                
                %     for i=1:length(KT)
                %
                %         y = prctile(KT{i}(mask_Fin),[1 99]);
                %
                %         KT{i} = EDTI_Library.E_DTI_outlier_smoothing(KT{i},mask_Fin,y);
                %
                %     end
                
            end
        end
        
        % From ExploreDTI: cast a 4D array into a 4D cell
        function B = E_DTI_DWI_mat2cell(A)
            
            if iscell(A)
                B = A;
            else
                for i=1:size(A,4)
                    B{i} = A(:,:,:,i);
                end
            end
        end
        
        % From ExploreDTI: cast a 4D cell into a 4D array
        function B = E_DTI_DWI_cell2mat(A)
            
            if iscell(A)
                B = repmat(A{1},[1 1 1 length(A)]);
                for i=1:length(A)
                    B(:,:,:,i) = A{i};
                end
            else
                B=A;
            end
        end
        
        % From ExploreDTI: local anisotropic smoothing
        function ret = E_DTI_smooth3_ani_voxel_size(data, filt, sz, arg, vox)
            
            if arg==0
                ret=data;
                return;
            end
            
            maskn = isnan(data);
            data(maskn) = 0;
            
            if length(sz)==1
                sz = [sz sz sz];
            end
            
            sz = sz(:)';
            
            padSize = (sz-1)/2;
            
            smooth = gaussian3(sz,arg,vox);
            
            ret=convn(EDTI_Library.padreplicate(data,padSize),smooth, 'valid');
            ret(maskn)=nan;
        end
        
        % From ExploreDTI: DTI/DKI fit - creates the .mat file - Adapted
        function E_DTI_model_fit(f_DWI,f_BM,fnam,Mask_par,perm,flip,fit_mode, dki_fit, dki_constraints, rekindle_kappa)
            global MRIToolkit;
            
            [DWI,VDims] = EDTI_Library.E_DTI_read_nifti_file(f_DWI);
            hdr = [];
            try
                hdr = load_untouch_header_only(f_DWI);
            catch
            end
            
            b = textread(f_BM);
            
            bval = sum(b(:,[1 4 6]),2);
            
            % Sort the acquisition on the fly
            bval = round(bval);
            [bval,IX] = sort(bval,'ascend');
            b = b(IX,:);
            DWI = DWI(:,:,:,IX);
            
            % The concept of b = 0 s/mm2 is a bit acquisition dependent.
            if(isfield(MRIToolkit,'min_bval_as_b0'))
                b0_indices = find(bval <= MRIToolkit.min_bval_as_b0);
            else
                b0_indices = find(bval <= 1);
                if(isempty(b0_indices))
                    min_bval = min(bval);
                    if(min_bval < 55) % bval tolerance set to 55, if we do not have any b = 0s/mm2
                       b0_indices = find(bval <= min_bval);
                    end
                end                
            end
%             if(~isempty(b0_indices) && isempty(intersect(b0_indices,1)))
%                warning('Found a b=0s/mm2, but this is not at the beginning. Please sort your data');               
%             end
%             order = diff(b0_indices);
%             discontinuity = find(order > 1,1);
%             if(~isempty(discontinuity))
%                 b0_indices = b0_indices(1:discontinuity);
%             end            
            NrB0 = length(b0_indices);
                        
            [bval, g] =  EDTI_Library.E_DTI_GetGradientsandBval_SC(b, NrB0);
            [g, b] = EDTI_Library.E_DTI_flip_g_b(g,b,perm,flip);
            
            bval = max(bval);
            
            par.clean_up_PIS = 1;
            par.ROBUST_option = 1;
            par.RE.rel_convergence = 1e-3;
            par.RE.max_iter = 20;
            if(exist('rekindle_kappa','var') > 0)
                par.RE.kappa = rekindle_kappa;
            else
                par.RE.kappa = 6;
            end
            par.DKI_constraints.do_it = 1;
            if(exist('dki_constraints','var') > 0)
                par.DKI_constraints.constr1 = dki_constraints;
            else
                par.DKI_constraints.constr1 = [-Inf Inf];
            end
            par.DKI_constraints.constr2 = [-Inf Inf];
            
            M = mean(DWI(:,:,:,NrB0+1:end),4);
            DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            
            if(isstruct(Mask_par))
                mask1 = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(DWI{1},Mask_par.tune_NDWI,Mask_par.mfs);
                mask2 = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(M,Mask_par.tune_DWI,Mask_par.mfs);
                
                mask = or(mask1,mask2);
            else
                mask = EDTI_Library.E_DTI_read_nifti_file(Mask_par);
                mask = mask > 0;
            end
            
            if(strcmpi(fit_mode,'ols'))
                par.TE = 1;
            elseif(strcmpi(fit_mode,'wls'))
                par.TE = 2;
            elseif(strcmpi(fit_mode,'rekindle'))
                par.TE = 4;
            elseif(strcmpi(fit_mode,'nls'))
                par.TE = 3;
            end
            
            g = [zeros(NrB0,3);g];
            
            tic;
            if dki_fit==0
                [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par,NrB0);
            elseif dki_fit==1
                [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par,NrB0,VDims);
            end
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
            
            g = g(NrB0+1:end,:);
            
            if ~isreal(FA)
                FA = real(FA);
            end
            if ~isreal(FEFA)
                FEFA = real(FEFA);
            end
            if ~isreal(FE)
                FE = real(FE);
            end
            if ~isreal(SE)
                SE = real(SE);
            end
            if ~isreal(eigval)
                eigval = real(eigval);
            end
            
            info = [];
            MDims = size(FA);
            
            try
                save(fnam,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
                    'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','par','Mask_par','hdr','-v7.3')
            catch
                save(fnam,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
                    'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','par','Mask_par','hdr')
            end
            
            if exist('KT','var') > 0
                save(fnam,'KT','-append')
            end
            
            toc;
        end
        
        % From ExploreDTI: A fast computation of the main diffusion directions
        function [V, D, PD] = E_DTI_eig_Lex_H(A)
            
            L = size(A,2);
            
            try
                V = repmat(double(0),[3 3 L]);
                D = repmat(double(0),[3 L]);
                A = double(A);
            catch
                V = repmat(single(0),[3 3 L]);
                D = repmat(single(0),[3 L]);
                A = single(A);
            end
            
            a1 = A(1,:);
            a2 = A(2,:);
            a3 = A(3,:);
            a4 = A(4,:);
            a5 = A(5,:);
            a6 = A(6,:);
            
            clear A;
            
            I1 = a1 + a4 + a6;
            I2 = a1.*a4 + a1.*a6 + a4.*a6;
            I2 = I2 - (a2.*a2 + a3.*a3 + a5.*a5);
            I3 = a1.*a4.*a6 + 2*a2.*a3.*a5;
            I3 = I3 - (a6.*a2.*a2 + a4.*a3.*a3 + a1.*a5.*a5);
            
            v = (I1.*I1)/9 - I2/3;
            s = (I1.*I1.*I1)/27 - (I1.*I2)/6 + I3/2;
            
            temp = ((s./v).*(1./sqrt(v)));
            temp(isnan(temp))=1;
            temp(temp>1)=1;
            temp(temp<-1)=-1;
            
            clear s;
            phi = acos(temp)/3;
            clear temp;
            
            D(1,:) = I1/3 + 2*sqrt(v).*cos(phi);
            D(2,:) = I1/3 - 2*sqrt(v).*cos(pi/3 + phi);
            D(3,:) = I1/3 - 2*sqrt(v).*cos(pi/3 - phi);
            
            
            PD = repmat(1>0,[4 L]);
            PD(1,:) = I3>0;
            clear I1 I2 I3 v phi;
            
            for i=1:3
                
                A = a1-D(i,:);
                B = a4-D(i,:);
                C = a6-D(i,:);
                
                V(1,i,:) = (a2.*a5 - B.*a3).*(a3.*a5 - C.*a2);
                V(2,i,:) = (a3.*a5 - C.*a2).*(a3.*a2 - A.*a5);
                V(3,i,:) = (a2.*a5 - B.*a3).*(a3.*a2 - A.*a5);
                
            end
            
            clear A B C;
            
            qv1 = V(:,1,:).*V(:,1,:);
            qv1 = sum(qv1,1);
            qv1 = sqrt(qv1);
            qv1 = repmat(qv1,[3 1 1]);
            qv1(qv1==0)=nan;
            V(:,1,:) = V(:,1,:)./qv1;
            clear qv1
            
            qv2 = V(:,2,:).*V(:,2,:);
            qv2 = sum(qv2,1);
            qv2 = sqrt(qv2);
            qv2 = repmat(qv2,[3 1 1]);
            qv2(qv2==0)=nan;
            V(:,2,:) = V(:,2,:)./qv2;
            clear qv2;
            
            qv3 = V(:,3,:).*V(:,3,:);
            qv3 = sum(qv3,1);
            qv3 = sqrt(qv3);
            qv3 = repmat(qv3,[3 1 1]);
            qv3(qv3==0)=nan;
            V(:,3,:) = V(:,3,:)./qv3;
            clear qv3;
            
            PD(2,:) = (a1.*a4 - a2.*a2)>=0;
            PD(3,:) = (a1.*a6 - a3.*a3)>=0;
            PD(4,:) = (a4.*a6 - a5.*a5)>=0;
            
            PD = all(PD,1);
            
            eps_a = 10^-4;
            
            % near degeneration
            dot_p1 = abs(squeeze(sum(V(:,1,:).*V(:,2,:),1)));
            dot_p2 = abs(squeeze(sum(V(:,2,:).*V(:,3,:),1)));
            dot_p3 = abs(squeeze(sum(V(:,1,:).*V(:,3,:),1)));
            dot_p = and(and(dot_p1(:)>eps_a,dot_p2(:)>eps_a),dot_p3(:)>eps_a);
            
            V(:,2,dot_p) = randn(3,1,sum(dot_p));
            V(:,2,dot_p) = (V(:,2,dot_p) - repmat(sum(V(:,2,dot_p).*V(:,1,dot_p),1),[3 1 1]).*V(:,1,dot_p));
            V(:,2,dot_p) = V(:,2,dot_p)./repmat(sqrt(sum(V(:,2,dot_p).*V(:,2,dot_p),1)),[3 1 1]);
            
            V(1,3,dot_p) = V(3,2,dot_p).*V(2,1,dot_p) - V(2,2,dot_p).*V(3,1,dot_p);
            V(2,3,dot_p) = V(1,2,dot_p).*V(3,1,dot_p) - V(3,2,dot_p).*V(1,1,dot_p);
            V(3,3,dot_p) = V(2,2,dot_p).*V(1,1,dot_p) - V(1,2,dot_p).*V(2,1,dot_p);
            
        end
        
        % From ExploreDTI: Compute the eigenvalues. Likely redundant.
        function [eigval, FE, SE, TE] = ComputeEigen(DT)
            eigval = zeros([3 size(DT,2)], 'double');
            FE = zeros([3 size(DT,2)], 'double');
            SE = zeros([3 size(DT,2)], 'double');
            TE = zeros([3 size(DT,2)], 'double');
            for i = 1:size(DT,2)
                [Vect,Val] = eig(reshape(DT(:,i),3,3));
                [eigval(:,i),idx] = sort(abs(diag(Val)),1,'descend');
                %     [eigval(:,i),idx] = sort(diag(Val),1,'descend');
                FE(:,i) = Vect(:,idx(1));
                SE(:,i) = Vect(:,idx(2));
                TE(:,i) = Vect(:,idx(3));
            end
        end
        
        % From ExploreDTI: Compute FA for a single set of Eigenvalues.
        function FA = ComputeFA(eigval)
            FA = sqrt(1/2).*sqrt((eigval(1,:)-eigval(2,:)).^2+(eigval(2,:)-eigval(3,:)).^2+(eigval(3,:)-eigval(1,:)).^2)./sqrt(eigval(1,:).^2+eigval(2,:).^2+eigval(3,:).^2);
        end
        
        % From ExploreDTI: FACT tractography from a seed point
        function tract = DTI_Track(param, seedPoint, seedDir)
            tract = cell(1,size(seedPoint,2));
            flist = 1:size(seedPoint,2);
            point = seedPoint;
            stepDir = seedDir;
            
            maxIt = param.maxLength/param.stepSize;
            it = 0;
            while it < maxIt
                it = it + 1;
                %     if param.pb; EDTI_Library.progressbar(size(point,2)); end;
                
                % update points
                point = point+param.stepSize.*stepDir;
                
                % get function of interest at current points
                DT = Interpolate(param.DT, point, param.VDims);
                
                % mask out all nan's
                mask0 = ~isnan(DT(1,:));
                DT = DT(:,mask0);
                point = point(:,mask0);
                stepDir = stepDir(:,mask0);
                flist = flist(mask0);
                
                % perform dti
                
                [eigval, newStepDir] = EDTI_Library.ComputeEigen(DT); clear DT;
                val = EDTI_Library.ComputeFA(eigval); clear eigval;
                
                % mask out small & NaN peaks
                mask2 = val > param.threshold;
                point = point(:,mask2);
                stepDir = stepDir(:,mask2);
                newStepDir = newStepDir(:,mask2);
                flist = flist(mask2);
                
                % mask out large angles
                a = EDTI_Library.angle2(stepDir,newStepDir);
                mask3 = a < param.maxAngle;
                point = point(:,mask3);
                stepDir = stepDir(:,mask3);
                newStepDir = newStepDir(:,mask3);
                flist = flist(mask3);
                
                % make sure we don't move back in the streamline
                flipsign = sign(sum(stepDir.*newStepDir,1));
                
                % update dir
                stepDir = flipsign([1 1 1],:).*newStepDir;
                
                % stop if we are out of points
                if isempty(point)
                    break
                end
                
                % add points to the tracts
                for i=1:length(flist)
                    tract{flist(i)}(it,:) = point(:,i);
                end
            end
        end
        
        % From ExploreDTI: The actual DTI-tractography trajectories calculation -
        % Adapted
        function Tracts = DTI_Tractography(param, DT, seedPoint, initSeedDir)
            
            
            DT = EDTI_Library.E_DTI_DT_cell2mat(DT);
            param.DT = single(DT); clear DT;
            
            % get function of interest at current points
            DT = Interpolate(param.DT, seedPoint, param.VDims);
            
            % mask out all nan's
            mask0 = ~isnan(DT(1,:));
            DT = DT(:,mask0);
            seedPoint = seedPoint(:,mask0);
            
            
            DT = DT([1 2 3 5 6 9],:);
            
            [V, eigval, PD] = EDTI_Library.E_DTI_eig_Lex_H(DT);
            
            clear DT PD;
            
            V = [squeeze(V(1,1,:))'; squeeze(V(2,1,:))'; squeeze(V(3,1,:))'];
            seedDir = V; clear V;
            
            
            val = EDTI_Library.ComputeFA(eigval); clear eigval;
            
            % enforce initSeedDir
            if nargin == 4
                flipsign = round(sum(initSeedDir.*seedDir,1));
                seedDir = repmat(flipsign,3,1).*seedDir;
            end
            
            % mask out all small or nan peaks
            mask2 = val > param.threshold;
            seedPoint = seedPoint(:,mask2);
            seedDir = seedDir(:,mask2);
            
            % if param.pb; EDTI_Library.progressbar('start', [size(seedPoint,2) 1], 'tracking in first direction'); end;
            
            Tracts1 = EDTI_Library.DTI_Track(param, seedPoint, seedDir);
            % if param.pb; EDTI_Library.progressbar('ready'); end;
            
            if nargin < 4
                %     if param.pb; EDTI_Library.progressbar('start', [size(seedPoint,2) 1], 'tracking in second direction'); end;
                Tracts2 = EDTI_Library.DTI_Track(param, seedPoint, -seedDir);
                %     if param.pb; EDTI_Library.progressbar('ready'); end;
            end
            
            % join Tracts
            if nargin < 4
                Tracts1 = cellfun(@flipud,Tracts1,'UniformOutput',false);
                Tracts = cell(size(Tracts1));
                for j = 1:size(Tracts1,2)
                    if ~isempty(Tracts1{j})
                        Tracts{j} = [Tracts1{j}; seedPoint(:,j)'];
                        if ~isempty(Tracts2{j})
                            Tracts{j} = [Tracts{j}; Tracts2{j}];
                        end
                    else
                        if ~isempty(Tracts2{j})
                            Tracts{j} = [seedPoint(:,j)'; Tracts2{j}];
                        end
                    end
                end
            else
                Tracts = cell(size(Tracts1));
                for j = 1:size(Tracts1,2)
                    if ~isempty(Tracts1{j})
                        Tracts{j} = [seedPoint(:,j)'; Tracts1{j};];
                    end
                end
            end
            
            % enforce length limitations
            % maska = cellfun('size',Tracts,1)*param.stepSize >= param.minLength;
            % maskb = cellfun('size',Tracts,1)*param.stepSize <= param.maxLength;
            % mask = maska & maskb;
            % Tracts = Tracts(mask);
            
            mask = cellfun('size',Tracts,1) > 1;
            Tracts = Tracts(mask);
            
        end
        
        % From ExploreDTI: Compute the linear, planar and spherical Westin's
        % measures
        function [CL, CP, CS] = E_DTI_Westin_measures(eigval)
            
            if ndims(eigval)==4
                
                L1 = eigval(:,:,:,1);
                L2 = eigval(:,:,:,2);
                L3 = eigval(:,:,:,3);
                
                clear eigval;
                
                mask = ~isnan(L1);
                
                CL = L1;
                CP = L1;
                CS = L1;
                
                CL(mask) = (L1(mask)-L2(mask))./L1(mask);
                CP(mask) = (L2(mask)-L3(mask))./L1(mask);
                CS(mask) = L3(mask)./L1(mask);
                
            else
                
                L1 = eigval(:,1);
                L2 = eigval(:,2);
                L3 = eigval(:,3);
                
                clear eigval;
                
                CL = (L1-L2)./L1;
                CP = (L2-L3)./L1;
                CS = L3./L1;
                
            end
        end
        
        % From ExploreDTI: Compute the diffusion measures given a set of fiber
        % tractography points
        function [t_FA, t_FE, t_Ang, t_GEO, t_L, t_MD] = E_DTI_diff_measures_vectorized(Tracts, VDims, DT)
            
            L = length(Tracts);
            
            Lts = zeros(L,1);
            
            for i=1:L
                Lts(i)=size(Tracts{i},1);
            end
            
            CLts = [0 ;cumsum(Lts)];
            
            Tmat = repmat(single(0),[sum(Lts) 3]);
            
            for i=1:L
                Tmat(CLts(i)+1:CLts(i+1),:) = Tracts{i};
            end
            
            Tmat = Tmat';
            
            DT = EDTI_Library.E_DTI_DT_cell2mat(DT);
            
            DT = Interpolate(single(DT), single(Tmat), single(VDims));
            
            DT = DT([1 2 3 5 6 9],:);
            
            [V, D, PD] = EDTI_Library.E_DTI_eig_Lex_H(DT);
            
            clear DT PD;
            
            V = [squeeze(V(1,1,:))'; squeeze(V(2,1,:))'; squeeze(V(3,1,:))'];
            
            dummy = EDTI_Library.FrAn_calc(D(1,:),D(2,:),D(3,:))*sqrt(3);
            
            dummy(dummy>sqrt(3)) = sqrt(3);
            dummy(dummy<0) = 0;
            
            t_FA = cell(1,L);
            
            for i=1:L
                t_FA{i} = dummy(1,CLts(i)+1:CLts(i+1))';
            end
            
            
            t_L = cell(1,L);
            for i=1:L
                t_L{i} = D(:,CLts(i)+1:CLts(i+1))';
            end
            
            dummy = sum(D,1)/3;
            
            t_MD = cell(1,L);
            
            for i=1:L
                t_MD{i} = dummy(1,CLts(i)+1:CLts(i+1))';
            end
            
            [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(D');
            
            clear D;
            
            t_GEO = cell(1,L);
            
            for i=1:L
                t_GEO{i} = [CL(CLts(i)+1:CLts(i+1)) CP(CLts(i)+1:CLts(i+1)) CS(CLts(i)+1:CLts(i+1))] ;
            end
            
            
            
            t_FE = EDTI_Library.E_DTI_Get_tract_dir(Tracts);
            t_Ang = cell(1,L);
            
            for i=1:L
                
                T = Tracts{i}(2:end,:)-Tracts{i}(1:end-1,:);
                T = T./repmat(sqrt(sum(T.*T,2)),[1 3]);
                
                t_FE{i} = abs(t_FE{i});
                
                t_Ang{i} = t_FA{i};
                
                t_Ang{i}(2:end-1)=(180/pi)*acos(abs(sum(T(1:end-1,:).*T(2:end,:),2)));
                t_Ang{i}(1) = (180/pi)*acos(abs(sum(T(1,:).*V(:,CLts(i)+1)',2)));
                t_Ang{i}(end) = (180/pi)*acos(abs(sum(T(end,:).*V(:,CLts(i+1))',2)));
                t_Ang{i} = real(t_Ang{i});
                
            end
        end
        
        % From ExploreDTI: helper for DTI-based tractography
        function Tract_dir = E_DTI_Get_tract_dir(Tracts)
            
            if length(Tracts)==0
                Tract_dir = cell(1,0);
            end
            
            for i=1:length(Tracts)
                
                Tract_dir{i} = Tracts{i};
                
                if size(Tracts{i},1)>2
                    dummy = Tract_dir{i}(2:end,:)-Tract_dir{i}(1:end-1,:);
                    Tract_dir{i}(2:end-1,:) = (dummy(1:end-1,:)+dummy(2:end,:))/2;
                    Tract_dir{i}(1,:) = dummy(1,:);
                    Tract_dir{i}(end,:) = dummy(end,:);
                else
                    dummy = Tract_dir{i}(2,:)-Tract_dir{i}(1,:);
                    Tract_dir{i}(1,:) = dummy;
                    Tract_dir{i}(2,:) = dummy;
                end
                
                
                NT = sqrt(sum(Tract_dir{i}.^2,2));
                
                Tract_dir{i} = Tract_dir{i}./[NT NT NT];
                
                
            end
            
        end
        
        % From ExploreDTI: DTI based fiber tractography - Adapted
        function WholeBrainTrackingDTI_fast(filename_in, filename_out, parameters)
            
            % parameters:
            %  SeedPointRes: [3 3 3]
            %             StepSize: 1
            %             FAThresh: 0.2000
            %          AngleThresh: 45
            %     FiberLengthRange: [50 500]
            
            
            try
                disp(['Loading ' filename_in '...']);
                warning off all
                load(filename_in,'VDims','FA','DT');
                warning on all
                disp(['Loading ' filename_in '... DONE!']);
                
                if ~exist('DT','var') || ~exist('FA','var') || ~exist('VDims','var')
                    disp('Error: incorrect or corrupted input file...');
                    return;
                end
                
            catch %ME
                disp('Error: incorrect or corrupted input file...');
                %     disp(['Error: ' ME.message]);
                return;
            end
            tic;
            
            DT = EDTI_Library.E_DTI_DT_cell2mat(DT);
            
            mask = ~isnan(FA); clear FA;
            
            disp('Calculating seed points...');
            seed_mask = mask; clear mask;
            
            S = size(seed_mask);
            x = single(VDims(1):VDims(1):VDims(1)*S(1));
            lx = length(x);
            y = single(VDims(2):VDims(2):VDims(2)*S(2));
            ly = length(y);
            z = single(VDims(3):VDims(3):VDims(3)*S(3));
            lz = length(z);
            
            X = repmat(x',[1 ly lz]);
            Y = repmat(y,[lx 1 lz]);
            p = zeros(1,1,lz);
            p(:,:,:) = z;
            Z = repmat(p,[lx ly 1]);
            clear x y z p lx ly lz;
            
            x_min = min(X(seed_mask));
            x_max = max(X(seed_mask));
            y_min = min(Y(seed_mask));
            y_max = max(Y(seed_mask));
            z_min = min(Z(seed_mask));
            z_max = max(Z(seed_mask));
            clear X Y Z
            
            xrange = x_min:parameters.SeedPointRes(1):x_max;
            yrange = y_min:parameters.SeedPointRes(2):y_max;
            zrange = z_min:parameters.SeedPointRes(3):z_max;
            
            seedPoint = zeros(3,size(xrange,2)*size(yrange,2)*size(zrange,2), 'single');
            i = 0;
            for x = xrange
                for y = yrange
                    for z = zrange
                        v = floor([x y z]./VDims);
                        if(seed_mask(v(1),v(2),v(3)))
                            i = i + 1;
                            seedPoint(:,i) = [x y z];
                        end
                    end
                end
            end
            seedPoint = seedPoint(:,1:i);
            disp('Calculating seed points...DONE!');
            clear seed_mask;
            % Tractography parameters
            param.VDims = single(VDims); % in mm
            param.stepSize = single(parameters.StepSize); % in mm
            param.maxAngle = single(parameters.AngleThresh); % in degrees
            param.threshold = single(parameters.FAThresh); % min FA value
            param.minLength = single(parameters.FiberLengthRange(1)); % min fiber length
            param.maxLength = single(parameters.FiberLengthRange(2)); % min fiber length
            param.pb = true;
            
            
            % Calculate trajectories
            disp('Calculating trajectories...');
            Tracts = EDTI_Library.DTI_Tractography(param, DT, seedPoint);
            disp('Calculating trajectories...DONE!');
            
            disp('Processing trajectories...');
            num_tracts = size(Tracts,2);
            FList = 1:num_tracts;
            
            [TractFA, TractFE, TractAng, TractGEO, TractLambdas, TractMD] =...
                EDTI_Library.E_DTI_diff_measures_vectorized(Tracts, param.VDims, DT);
            
            TractL = cell(1,num_tracts);
            
            for i = 1:num_tracts
                TractL{i} = single(size(Tracts{i},1));
            end
            
            
            TractMask = repmat(uint16(0),S);
            
            for i = FList
                IND = unique(sub2ind(S,...
                    round(double(Tracts{i}(:,1))./VDims(1)),...
                    round(double(Tracts{i}(:,2))./VDims(2)),...
                    round(double(Tracts{i}(:,3))./VDims(3))));
                TractMask(IND) = uint16(double(TractMask(IND)) + 1);
            end
            
            clear DT;
            disp('Processing trajectories...DONE!');
            
            % Save everything
            disp('Saving trajectories...')
            save(filename_out,'FList'); clear FList;
            save(filename_out,'Tracts','-append'); clear Tracts;
            save(filename_out,'TractL','-append'); clear TractL;
            save(filename_out,'TractFE','-append'); clear TractFE;
            save(filename_out,'TractFA','-append'); clear TractFA;
            save(filename_out,'TractAng','-append'); clear TractAng;
            save(filename_out,'TractGEO','-append'); clear TractGEO;
            save(filename_out,'TractMD','-append'); clear TractMD;
            save(filename_out,'TractLambdas','-append'); clear TractLambdas;
            save(filename_out,'TractMask','VDims','-append'); clear TractMask;
            disp('Saving trajectories...DONE!')
            t=toc;
            
            m=t/60;
            
            disp(['Total computation time was ' num2str(m) ' min.'])
            
        end
        
        % From ExploreDTI: Average dwi signals within an FA thresholded mask
        function dwi = E_DTI_get_dwi_4_resp_func_2(DWI, FA_T, FA)
            
            if ~iscell(DWI)
                DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            end
            
            if FA_T~=0
                fa_mask = FA>=FA_T*sqrt(3);
            elseif FA_T==0
                fa_mask = ~isnan(FA);
            end
            
            dwi = zeros(length(DWI),sum(fa_mask(:)),'single');
            
            for i=1:length(DWI)
                dwi(i,:) = single(DWI{i}(fa_mask(:))');
            end
        end
        
        % From ExploreDTI: creates an RF given some initial parameters
        function shcoef = E_DTI_create_initial_RF_RC(lmax,bval,FA,trD)
            
            file=load('icosahedron5.mat');
            g=file.Expression1;
            b=repmat(bval,[length(g) 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
            [lambda1,lambda2,lambda3] = EDTI_Library.Calculate_lambdas(FA,trD,1);
            Drot=(EDTI_Library.Rezgamma(0)*EDTI_Library.Reybeta(0))*([lambda2 0 0;0 lambda2 0;0 0 lambda1])*((EDTI_Library.Rezgamma(0)*EDTI_Library.Reybeta(0))');
            s2=exp(-b*[Drot(1,1) Drot(1,2) Drot(1,3) Drot(2,2) Drot(2,3) Drot(3,3)]');
            
            sh = SH(lmax,g);
            shcoef = sh.coef(s2);
            j = 0;
            for l_ = 0:2:lmax
                for m = -l_:l_
                    j = j + 1;
                    if m ~= 0
                        shcoef(j,:) = 0;
                    end
                end
            end
            shcoef=shcoef/shcoef(1);
        end

        % Original: creates an RF given a single fiber signal
        function shcoef = E_DTI_create_initial_RF_RC_from_signal(lmax,s2)

            file=load('icosahedron5.mat');
            g=file.Expression1;
                        
            sh = SH(lmax,g);
            shcoef = sh.coef(s2);
            j = 0;
            for l_ = 0:2:lmax
                for m = -l_:l_
                    j = j + 1;
                    if m ~= 0
                        shcoef(j,:) = 0;
                    end
                end
            end
            shcoef=shcoef/shcoef(1);
        end
        

        % Adapted: simulate a single with the DTI model + isotropic
        % kurtosis
        function s2 = SimulateSignalDTIIsotropicK(bvals,g,eigenvalues,MK,angles)
            
            b=repmat(bvals,[1 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
            lambda1 = eigenvalues(1);
            lambda2 = eigenvalues(2);
            lambda3 = eigenvalues(3);
            MD = mean(eigenvalues);
            
            Drot=(EDTI_Library.Rezgamma(angles(1))*EDTI_Library.Reybeta(angles(2)))*([lambda3 0 0;0 lambda2 0;0 0 lambda1])*((EDTI_Library.Rezgamma(angles(1))*EDTI_Library.Reybeta(angles(2)))');
            s2=exp(-b*[Drot(1,1) Drot(1,2) Drot(1,3) Drot(2,2) Drot(2,3) Drot(3,3)]'+1/6*bvals.*bvals*MD*MD*MK);
        end        
        
        % From ExploreDTI: Helper of E_DTI_create_initial_RF_RC
        function [lambda1,lambda2,lambda3]=Calculate_lambdas(FA,trD,fractionlambda)
            
            if fractionlambda==1
                lambda1=(trD/3)*(1+2*FA/(3-2*FA^2)^(1/2));
                lambda2=(trD/3)*(1-FA/(3-2*FA^2)^(1/2));
                lambda3=lambda2;
            elseif fractionlambda<1&&fractionlambda>0
                lambda1=(1+fractionlambda^2+sqrt((1-2*fractionlambda-2*FA^2+fractionlambda^2*(1-2*FA^2))/(-3+2*FA^2))+fractionlambda*sqrt((1-2*fractionlambda-2*FA^2+fractionlambda^2*(1-2*FA^2))/(-3+2*FA^2)))*trD/(2*(1+fractionlambda+fractionlambda^2));
                lambda2=(1+fractionlambda-sqrt((1-2*fractionlambda-2*FA^2+fractionlambda^2*(1-2*FA^2))/(-3+2*FA^2)))*trD/(2*(1+fractionlambda+fractionlambda^2));
                lambda3=fractionlambda*(1+fractionlambda-sqrt((1-2*fractionlambda-2*FA^2+fractionlambda^2*(1-2*FA^2))/(-3+2*FA^2)))*trD/(2*(1+fractionlambda+fractionlambda^2));
            end
        end
        
        % From ExploreDTI: a rotation matrix around the Z axis
        function Rz=Rezgamma(gamma)
            Rz=[cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0;0 0 1];
        end
        
        % From ExploreDTI: a rotation matrix around the Y axis
        function Ry=Reybeta(beta)
            Ry=[cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)];
        end
        
        % From ExploreDTI: CSD with tensor based response function - Adapted
        function D = E_DTI_Get_HARDI_CSD_FOD(fin, sim_rf, Lmax, sim_adc, sim_fa, filename_out, save_sh)
            
            fn = [fin(1:end-4) '_CSD_FOD.nii'];
            
            if exist(fn,'file')==2
                [D, VDims] = EDTI_Library.E_DTI_read_nifti_file(fn);
                D = single(D);
                D(D==0)=nan;
                return;
            end
            
            disp('CSD with tensor based RF calibration (not recommended!)...')
            
            warning off all
            load(fin,'DWI','NrB0','b','g','FA','bval','MDims','DT','VDims')
            warning on all
            
            if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('b','var') || ...
                    ~exist('FA','var') || ~exist('g','var') || ~exist('bval','var') || ~exist('DT','var')
                disp(['Format of ''' fin ''' not correct, skipping data!'])
                D = [];
                return;
            end
            clear DT;
            
            mask = repmat(true,MDims);
            
            FA(~mask)=nan;
            
            bvals=round(sum(b(:,[1 4 6]),2)/100)*100;
            ubvals = unique(bvals);
            
            % bvals = repmat(bval,[length(DWI)-NrB0 1]);
            grad4 = [[zeros(NrB0,3);g] bvals];
            
            if length(ubvals)>2
                disp(['Selecting a subset of the data: ' num2str(ubvals(1)) 's/mm2 and ' num2str(ubvals(end)) 's/mm2']);
                IX = find(bvals <= ubvals(1) | bvals >= ubvals(end));
                DWI = DWI(IX);
                grad4 = grad4(IX,:);
                b = b(IX,:);
                bvals = bvals(IX);
            end
            
            DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            
            B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, ~isnan(FA), NrB0);
            DWI(1:NrB0)=B0s;
            clear B0s;
            
            if(sim_rf == 1)
                the_dirs = load('dir300.txt');
                s_grad4 = [the_dirs ...
                    repmat(bval,...
                    [size(the_dirs,1) 1])];
                s_grad4 = [0 0 0 0; s_grad4];
                
                dti = DTI(s_grad4);
                dwi = dti.sim_dwi(sim_fa,...
                    sim_adc,0,0,...
                    bval);
                r_sh = SD.response(dwi, s_grad4,...
                    bval, 0, Lmax);
            else
                
                dwi = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(DWI,0.7,FA);
                dwi(dwi<=0)=1;
                
                if numel(dwi)==0
                    D = [];
                    disp('Error calculating response function (DWIs)..., skipping data!')
                    return;
                end
                
                r_sh = SD.response(dwi, grad4,...
                    bval, 0.7, Lmax);
                
                if any(isnan(r_sh))
                    D = [];
                    disp('Error calculating response function..., skipping data!')
                    return;
                end
            end
            
            try
                sh = SH(Lmax,grad4(bvals>=ubvals(end),1:3));
                csd = CSD(r_sh,Lmax);
            catch
                D = [];
                disp('Error calculating CSD FOD..., skipping data!')
                return;
            end
            
            dwi = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(DWI,0,FA);
            dwi = dwi(NrB0+1:end,:);
            dwi(dwi<=0)=1;
            
            dwi_sh = sh.coef(dwi);
            
            if(exist('save_sh','var') > 0 && save_sh == 1)
                tmp_dwi_sh = unvec(dwi_sh,mask);
                EDTI_Library.E_DTI_write_nifti_file(tmp_dwi_sh,VDims,[fin(1:end-4) '_SHcoeffs.nii']);
                clear tmp_dwi_sh;
            end
            
            try
                D = zeros(size(dwi_sh),'single');
                %     E_DTI_open_matlabpool;
                parfor i=1:size(D,2)
                    D(:,i) = csd.deconv(dwi_sh(:,i));
                end
                %     E_DTI_close_matlabpool;
            catch
                D = csd.deconv(dwi_sh);
            end
            
            clear dwi_sh;
            
            D = unvec(D,~isnan(FA));
            D = single(D);
            
            
            DD = single(D);
            DD(isnan(DD))=0;
            
            % fn = [fin(1:end-4) '_CSD_FOD.nii'];
            
            if(exist('filename_out','var') > 0 && ~isempty(filename_out))
                fn = filename_out;
            else
                fn = [fin(1:end-4) '_CSD_FOD.nii'];
            end
            
            EDTI_Library.E_DTI_write_nifti_file(DD,VDims,fn);
        end
        
        % From ExploreDTI: CSD with recursively calibrated response function (Tax
        % et al.) - adapted
        function D = E_DTI_Get_HARDI_CSD_FOD_RC(fin,Lmax,rc_mask_file,file_out,save_sh)
            
            fn = file_out;%[fin(1:end-4) '_CSD_FOD.nii'];
            
            if exist(fn,'file')==2
                [D, VDims] = EDTI_Library.E_DTI_read_nifti_file(fn);
                D = single(D);
                D(D==0)=nan;
                return;
            end
            
            disp('CSD with recursive RF calibration...')
            
            
            if(exist('rc_mask_file','var') && ~isempty(rc_mask_file))
                disp(['Constraining within ' rc_mask_file]);
                rc_mask = EDTI_Library.E_DTI_read_nifti_file(rc_mask_file);
            end
            
            try
                pctRunOnAll warning off all
            catch
                warning off all
            end
            
            t = 0.01;
            lmax = Lmax;
            it=10;
            
            suf = 5;
            
            load(fin,'DWI','VDims','NrB0','b','g','FA','bval','MDims')
            
            if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('b','var') || ...
                    ~exist('FA','var') || ~exist('g','var') || ~exist('bval','var')
                disp(['Format of ''' fin ''' not correct, skipping data!'])
                D = [];
                return;
            end
            
            DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            
            B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, ~isnan(FA), NrB0);
            DWI(1:NrB0)=B0s;
            clear B0s;
            
            mask = ~isnan(FA);
            bvals=round(sum(b(:,[1 4 6]),2)/100)*100;
            ubvals = unique(bvals);
            
            grad=[[zeros(NrB0,3);g],bvals];
            
            if length(ubvals)>2
                disp(['Selecting a subset of the data: ' num2str(ubvals(1)) 's/mm2 and ' num2str(ubvals(end)) 's/mm2']);
                IX = find(bvals <= ubvals(1) | bvals >= ubvals(end));
                DWI = DWI(IX);
                grad = grad(IX,:);
                b = b(IX,:);
                bvals = bvals(IX);
            end
            
            dwi=vec(EDTI_Library.E_DTI_DWI_cell2mat(DWI),mask); clear DWI;
            dwi = single(dwi);
            
            M=FA/sqrt(3);
            M=vec(M,mask);
            
            % l=true(size(b,1),1);
            % dwi = double(dwi(l,:));
            % grad = grad(l,:);
            
            % Select voxels
            sf_mask=(double(M)>0.01); clear M;
            
            fim = find(sf_mask==1);
            
            if suf>length(fim)/500
                suf = round(length(fim)/500);
                if suf==0
                    suf=1;
                end
            end
            
            fim = fim(1:suf:end);
            sf_mask(:)=0;
            sf_mask(fim)=1;
            
            if(exist('rc_mask_file','var') && ~isempty(rc_mask_file))
                fprintf('Before sf had %d candidates %s',length(find(sf_mask>0)),newline);
                sf_mask = sf_mask & vec(rc_mask>0,mask);
                fprintf('After sf has %d candidates %s',length(find(sf_mask>0)),newline);
            end
            
            dwi = dwi(:,sf_mask);
            
            % select single shell without b0
            dwi = dwi(NrB0+1:end,:);
            grad = grad(NrB0+1:end,1:3);
            
            % dwi in SH
            sh = SH(lmax,grad);
            
            dwi_sh = sh.coef(dwi);
%             dwi_sh = sh.reg_coef(dwi,0.001);
            
            % axial symmetric RF in z direction for every voxel, FA tensor = 0.05,
            % werkt voor lmax = 6 en 8
            r_sh = zeros(SH.lmax2n(lmax),size(dwi,2));
            
            shcoef = EDTI_Library.E_DTI_create_initial_RF_RC(lmax,bval,0.05,2.1*10^(-3));
            
            for i=1:length(shcoef(:))
                
                r_sh(i,:) = shcoef(i,1);
                
            end
            
            SHPrecomp.init(lmax); % always perform this initialization first
            
            I=cell(1);
            for i=1:it
                %     disp(num2str(i))
                % total RF
                I{i}.nvox=size(r_sh,2); %amount of voxels
                r_shtot = mean(r_sh,2);
                r_shtot_sd = std(r_sh,0,2);
                I{i}.r_shtot=r_shtot;
                I{i}.r_shtot_sd = r_shtot_sd;
                
                % create CSD object
                csd = CSD(r_shtot,lmax); clear r_shtot
                
                % perform constrained deconvolution on SH coefficients
                csd_fod_sh = csd.deconv(dwi_sh);
                
                % peak extraction
                [dirs_, vals_] = SHPrecomp.all_peaks(csd_fod_sh, 0, 2);
                
                vals = zeros(2,length(vals_));
                
                for k=1:length(vals_)
                    if isempty(vals_{k})
                        vals_{k}(1) = nan;
                    end
                    vals(1,k) = vals_{k}(1);
                    if length(vals_{k})>1
                        vals(2,k) = vals_{k}(2);
                    else
                        vals(2,k) = nan;
                    end
                end
                
                dirs = zeros(3,length(vals_));
                
                for k=1:length(vals_)
                    if isempty(dirs_{k})
                        dirs_{k} = [nan nan nan]';
                    end
                    dirs(1,k) = dirs_{k}(1,1);
                    dirs(2,k) = dirs_{k}(2,1);
                    dirs(3,k) = dirs_{k}(3,1);
                end
                
                I{i}.vals=vals;
                peak_mask = ((vals(2,:)./vals(1,:)<t | isnan(vals(2,:))) & ~isnan(vals(1,:)));
                I{i}.peak_mask = peak_mask;
                dwi_sh = dwi_sh(:,peak_mask);
                r_sh=zeros(SH.lmax2n(lmax),size(dwi_sh,2));
                dirs=dirs(:,peak_mask);
                dwi=dwi(:,peak_mask);
                
                for j = 1:size(dwi_sh,2)
                    fe=dirs(1:3,j);
                    nullsp=null(fe');
                    se=nullsp(:,1);
                    te=nullsp(:,2);
                    rot_grad = zeros(size(grad,1),3);
                    rot = [te se fe];
                    for k = 1:size(grad,1)
                        rot_grad(k,:) = grad(k,:)*rot;
                    end
                    r_sh(:,j) = (SH.eval(lmax,c2s(rot_grad))\dwi(:,j));
                end
                
                j = 0;
                for l_ = 0:2:lmax
                    for m = -l_:l_
                        j = j + 1;
                        if m ~= 0
                            r_sh(j,:) = 0;
                        end
                    end
                end
                
                if i>1
                    change=abs(I{i}.r_shtot-I{i-1}.r_shtot)./I{i-1}.r_shtot;
                    if  all(change(~isnan(change))<0.01)
                        break
                    end
                end
                
                
            end
            
            % Plot mean PR and nr of voxels over iteration
            % nvox=[];
            % ratiovals=[];
            r_shtot=[];
            r_shtot_sd = [];
            for i=1:size(I,2)
                %     nvox=[nvox,I{i}.nvox];
                %     ratiovals=[ratiovals,mean([I{i}.vals(2,~isnan(I{i}.vals(2,:)))/I{i}.vals(1,~isnan(I{i}.vals(2,:))),zeros(1,length(I{i}.vals(2,isnan(I{i}.vals(2,:)))))])];
                r_shtot=[r_shtot,I{i}.r_shtot];
                r_shtot_sd=[r_shtot_sd,I{i}.r_shtot_sd];
            end
            
            % Save the Response function per iteration (ADL)
            save([fin(1:end-4) '_r_sh.mat'],'r_shtot','r_shtot_sd','-v7.3');
            
            load(fin,'DWI','VDims','NrB0','b','g','FA','bval','MDims')
            bvals=round(sum(b(:,[1 4 6]),2)/10)*10;
            
            if length(ubvals)>2
                disp(['Selecting a subset of the data: ' num2str(ubvals(1)) 's/mm2 and ' num2str(ubvals(end)) 's/mm2']);
                IX = find(bvals <= ubvals(1) | bvals >= ubvals(end));
                DWI = DWI(IX);
                %     grad = grad(IX,:);
                b = b(IX,:);
                bvals = bvals(IX);
            end
            %
            % % bvals=sum(b(:,[1 4 6]),2);
            % bvals=round(sum(b(:,[1 4 6]),2)/100)*100;
            % grad=[[zeros(NrB0,3);g],bvals];
            dwi=vec(EDTI_Library.E_DTI_DWI_cell2mat(DWI),mask); clear DWI;
            dwi = single(dwi);
            %
            % % Do CSD for iteration(s) iter, mostly last iteration
            iter=size(I,2);
            % % select single shell
            dwi = dwi(NrB0+1:end,:);
            % b = b(NrB0+1:end,:);
            % grad = grad(NrB0+1:end,1:3);
            % % create SH object
            % %%sh = SH(grad, lmax);
            % sh = SH(lmax, grad);
            %
            % % estimate SH coefficients
            dwi_sh = sh.coef(dwi);
            
            if(exist('save_sh','var') > 0 && save_sh == 1)
                tmp_dwi_sh = unvec(dwi_sh,mask);
                EDTI_Library.E_DTI_write_nifti_file(tmp_dwi_sh,VDims,[fin(1:end-4) '_SHcoeffs.nii']);
                clear tmp_dwi_sh;
            end
            
            for i=iter
                %     disp(num2str(i))
                r_sh = r_shtot(:,i);
                
                if any(isnan(r_sh))
                    D=[];
                    disp('Error during CSD. The estimated RF contains NaN values.');
                    return;
                end
                
                % create CSD object
                csd = CSD(r_sh,lmax);
                % perform constrained deconvolution on SH coefficients
                csd_fod_sh = csd.deconv(dwi_sh);
                D = unvec(csd_fod_sh,mask);
                %     save([directory,tag,'_CSD_it',num2str(i),'_fod.mat'],'csd_fod_sh');
                %     EDTI_Library.E_DTI_write_nifti_file(csd_fod_sh,VDims,[directory,tag,'_CSD_it',num2str(i),'_fod.nii'])
            end
            
            % For quality check, compute the residuals (ADL)
            residuals = mean(abs(sh.amp(csd.conv(csd_fod_sh))-dwi));
            residuals = unvec(residuals,mask);
            save([fin(1:end-4) '_csd_residuals.mat'],'residuals');
            
            DD = single(D);
            DD(isnan(DD))=0;
            
            if(exist('file_out','var') > 0 && ~isempty(file_out))
                fn = file_out;
            else
                if(exist('rc_mask_file','var') > 0 && ~isempty(rc_mask_file))
                    fn = [fin(1:end-4) '_' rc_mask_file];
                else
                    fn = [fin(1:end-4) '_CSD_FOD.nii'];
                end
            end
            
            EDTI_Library.E_DTI_write_nifti_file(DD,VDims,fn);
            
            try
                pctRunOnAll warning on all
            catch
                warning on all
            end
        end
        
        % From ExploreDTI: Multi Shell CSD (MSCSD, Jeurissen et al.) - adapted
        function D = E_DTI_Get_HARDI_CSD_FOD_RC_MuSh(file_in,lmax,rc_mask_file,file_out,fast_pve_seg_file)
            
            try
                pctRunOnAll warning off all;
            catch
                warning off all;
            end
            
            par.wm_fa = 0.7;
            par.gm_fa = 0.2;
            par.cs_fa = 0.2;
            
            par.tau = 0.1;
            par.nitts  = 50;
            par.lmax = lmax;
            par.lambda = 0.1;
            
            load(file_in, 'DWI', 'NrB0', 'b', 'g', 'FA', 'eigval','VDims');
            
            bvals = sum(b(:,[1 4 6]),2);
            IX = find(bvals <= min(bvals) | (bvals > 200 & bvals < 1300));
            
            mask = ~isnan(FA);
            
            par.clean_up_PIS = 1;
            par.RE.rel_convergence = 1e-3;
            par.RE.max_iter = 20;
            par.RE.kappa = 6;
            par.ROBUST_option = 1;
            
            par.TE = 1;
            
            DT = EDTI_Library.E_DTI_Get_DT_from_DWI_b_mask(DWI(IX),b(IX,:),mask,par,NrB0);
            
            FA = FA./sqrt(3);
            FA(isnan(FA))=0; FA_thresh = 0.00001;
            
            dwi = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(DWI, FA_thresh, FA);
            m = mean(dwi(1:NrB0,:),1);
            dwi = dwi ./ repmat(m/100, [size(dwi, 1) 1]);
            DWI = unvec(dwi, FA>FA_thresh);
            
            [f_pn,nm] = fileparts(file_in); load(file_in, 'VDims');
            if(isempty(f_pn))
                f_pn = pwd;
            end
            
            if(exist('fast_pve_seg_file','var') && ~isempty(fast_pve_seg_file))
                disp('Using provided T1 segmentation to extract the response functions');
                A = EDTI_Library.E_DTI_read_nifti_file(fast_pve_seg_file);
                wm_mask = A == 3;
                gm_mask = A == 2;
                cs_mask = A == 1;
            else
                disp('Deriving segmentation from dMRI - heuristic (use T1 if available)');
                [~, FA_c, ~, ~, eigval_c] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
                MD_c = sum(eigval_c,4)/3;
                wm_mask = FA_c > 0.2 & MD_c < 0.001;
                gm_mask = FA_c < 0.2 & MD_c < 0.001;
                cs_mask = MD_c > 0.001 & FA_c < 0.2;
            end
            
            EDTI_Library.E_DTI_write_nifti_file(single(wm_mask), VDims, [f_pn filesep nm '_WM_mask.nii']);
            EDTI_Library.E_DTI_write_nifti_file(single(gm_mask), VDims, [f_pn filesep nm '_GM_mask.nii']);
            EDTI_Library.E_DTI_write_nifti_file(single(cs_mask), VDims, [f_pn filesep nm '_CSF_mask.nii']);
            
            [D, max_kernel] = EDTI_Library.E_DTI_MS_CSD_deconvolution_MuSh_mask(file_in, par, wm_mask, gm_mask, cs_mask, FA_thresh, rc_mask_file);
            
            lmax = EDTI_Library.E_DTI_n2lmax(max_kernel);
            
            gm_volume = D(SH.lmax2n(lmax)+1, :);
            cs_volume = D(SH.lmax2n(lmax)+2, :);
            
            GM_volume = unvec(gm_volume, FA>FA_thresh);
            CS_volume = unvec(cs_volume, FA>FA_thresh);
            
            D = D(1:SH.lmax2n(lmax), :);
            D = unvec(D, FA>FA_thresh);
            
            WM_volume = D(:,:,:,1);
            
            TOT = WM_volume + GM_volume + CS_volume;
            WM_volume = WM_volume ./ TOT;
            GM_volume = GM_volume ./ TOT;
            CS_volume = CS_volume ./ TOT;
            
            WM_volume(~isfinite(WM_volume)) = 0;
            GM_volume(~isfinite(GM_volume)) = 0;
            CS_volume(~isfinite(CS_volume)) = 0;
            
            D = single(D);
            D(isnan(D))=0;
            
            % disp('saving state for debug');
            EDTI_Library.E_DTI_write_nifti_file(D, VDims, file_out);
            EDTI_Library.E_DTI_write_nifti_file(WM_volume, VDims, [file_out(1:end-4) '_WM_fraction.nii']);
            EDTI_Library.E_DTI_write_nifti_file(GM_volume, VDims, [file_out(1:end-4) '_GM_fraction.nii']);
            EDTI_Library.E_DTI_write_nifti_file(CS_volume, VDims, [file_out(1:end-4) '_CSF_fraction.nii']);
            
            try
                pctRunOnAll warning on all;
            catch
                warning on all;
            end
            
        end
        
        % From ExploreDTI: Perform the actual multi-shell deconvolution
        function [ CSD_FOD, max_kernel] = E_DTI_MS_CSD_deconvolution_MuSh_mask(file, par, wm_mask, gm_mask, cs_mask, FA_thresh,sf_mask)
            
            load(file, 'DWI', 'g', 'FA', 'NrB0','b') % 'bvals',
            
            if(exist('sf_mask','var') && ~isempty(sf_mask))
                if(ischar(sf_mask))
                    sf_mask = EDTI_Library.E_DTI_read_nifti_file(sf_mask);
                end
                wm_mask = wm_mask & sf_mask;
            end
            
            % vectorise the data
            DWI = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(DWI, FA_thresh, FA);
            wm_mask = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(wm_mask, FA_thresh, FA);
            gm_mask = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(gm_mask, FA_thresh, FA);
            cs_mask = EDTI_Library.E_DTI_get_dwi_4_resp_func_2(cs_mask, FA_thresh, FA);
            wm_mask = logical(wm_mask);
            gm_mask = logical(gm_mask);
            cs_mask = logical(cs_mask);
            
            % by f.guo
            bvals = sum(b(:,[1 4 6]),2);
%             bvals = round(bvals/10) *10;
            bvals = round(bvals/100) *100;
            
            % determine the number of unique shells and order by increasing b-value
            bvalues = unique(bvals);
            bvalues = sort(bvalues, 'ascend');
            nshells = length(bvalues);
            
            % break the data down into individual shells
            % and, where b>0, grab the g-matricies for those shells
            shells = cell(1, nshells);
            g_mats = cell(1, nshells);
            btrunk = bvals(NrB0+1:size(DWI,1));
            for x=1:nshells
                shells{x} = DWI(bvals==bvalues(x), :);
                if x==1
                    g_mats{x} = zeros(NrB0, 3);
                else
                    g_mats{x} = g(btrunk==bvalues(x), :);
                end
            end
            
            % determine the maximum lmax supported by each shell
            lmax = zeros(1, nshells);
            for x=1:nshells
                lmax(x) = EDTI_Library.E_DTI_n2lmax(size(shells{x}, 1));
            end
            
            lmax = floor(lmax/2)*2;
            
            % if the user manually specified a value for each shell, check those values
            % are supported
            if numel(par.lmax)==length(shells)
                if any(par.lmax > lmax)
                    warning(['One or more user supplied shell-wise lmax is unsupported',...
                        ' by the number of datapoints in the corresponding shell.']);
                    warning(['Changing : ' num2str(par.lmax)]);
                    par.lmax(par.lmax > lmax) = lmax(par.lmax > lmax);
                    warning(['To       : ' num2str(par.lmax)]);
                end
                if par.lmax(1) ~= 0
                    warning('User specified lmax for b0 shell is non-zero. Making alteration');
                    par.lmax(1) = 0;
                end
                lmax = par.lmax;
            else
                % if the user specified a single lmax, use that as the maximum
                if numel(par.lmax)==1 && par.lmax > 0
                    lmax(lmax>par.lmax)=par.lmax;
                end
                
                % set the first lmax (the one used on the b0 shell) to zero
                lmax(1) = 0;
            end
            
            % now calculate response functions for each shell
            r_sh = cell(nshells, 3);
            conv = cell(nshells, 3);
            
            sel_voxels = zeros(size(shells{1},2),nshells);
            
            for x=1:nshells
                %     disp(['Fitting response functions to shell ' num2str(x) ' of ' num2str(nshells)]);
                wm_vox = shells{x}(:, wm_mask);
                gm_vox = shells{x}(:, gm_mask);
                cs_vox = shells{x}(:, cs_mask);
                
                % slightly different proceduces based on whether this is the b0 shell
                % or the shells with b>0
                if x==1 || lmax(x) == 0
                    % where the b0 is involved, just fit an isotropic response
                    r_sh{x, 1} = EDTI_Library.E_DTI_MS_b0_response_MuSh(wm_vox, size(wm_vox,1));
                    r_sh{x, 2} = EDTI_Library.E_DTI_MS_b0_response_MuSh(gm_vox, size(wm_vox,1));
                    r_sh{x, 3} = EDTI_Library.E_DTI_MS_b0_response_MuSh(cs_vox, size(wm_vox,1));
                    
                else
                    % otherwise, need to do some g-matrix reorientation
                    
                    % set up the grad matricies for tensor based reorientation
                    g = g_mats{x};
                    bvals = repmat(bvalues(x), [size(wm_vox, 1) 1]);
                    grad4 = [zeros(NrB0, 4);[g bvals]];
                    
                    % need to stack on the b0's for tensor based reorientation
                    wm_vox = [shells{1}(:, wm_mask); wm_vox];
                    gm_vox = [shells{1}(:, gm_mask); gm_vox];
                    cs_vox = [shells{1}(:, cs_mask); cs_vox];
                    
                    % chop out a representative sample
                    
                    [wm_vox,wm_vox_subsamp] = EDTI_Library.UniformMaskSubsamp(wm_vox,5000);
                    gm_vox = EDTI_Library.UniformMaskSubsamp(gm_vox,5000);
                    cs_vox = EDTI_Library.UniformMaskSubsamp(cs_vox,5000);
                    
                    %         wm_vox = wm_vox(:, 1:min(size(wm_vox, 2), 5000));
                    %         gm_vox = gm_vox(:, 1:min(size(gm_vox, 2), 5000));
                    %         cs_vox = cs_vox(:, 1:min(size(cs_vox, 2), 5000));
                    
                    % and fit the response functions
                    %       r_sh{x, 1} = SD.response(wm_vox, grad4, bvalues(x), 0.7, lmax(x));
                    
                    
                    [r_sh{x, 1},sf_voxels] = EDTI_Library.E_DTI_HARDI_CSD_FOD_RC_MuSh_mask(wm_vox(NrB0+1:end,:),grad4(NrB0+1:end,:),bvalues(x),lmax(x));
                    r_sh{x, 2} = SD.response(gm_vox, grad4, bvalues(x), 0, 0);
                    r_sh{x, 3} = SD.response(cs_vox, grad4, bvalues(x), 0, 0);
                    
                    idx_conv = find(wm_mask);
                    sel_voxels(idx_conv(wm_vox_subsamp),x) = sf_voxels;
                    
                end
                
                % then use those responses to build convolution kernels
                conv{x, 1} = EDTI_Library.E_DTI_calc_convmat_MuSh(r_sh{x, 1}, lmax(x));
                conv{x, 2} = EDTI_Library.E_DTI_calc_convmat_MuSh(r_sh{x, 2}, 0);
                conv{x, 3} = EDTI_Library.E_DTI_calc_convmat_MuSh(r_sh{x, 3}, 0);
                
            end
            
            sel_voxels_vol = unvec(sel_voxels',FA>FA_thresh);
            save('MSCSD_SF_voxels','sel_voxels_vol');
            save('MSCSD_sh_RCresponse','r_sh')
            
            % then build a combined convolution kernel
            fconv = [];
            max_kernel  = SH.lmax2n(max(lmax)); % the largest convolution kernel size
            for x=1:nshells
                
                % extract the convolution kernels
                conv_wm = conv{x, 1};
                conv_gm = conv{x, 2};
                conv_cs = conv{x, 3};
                
                % combine them into a convolution block for the shell in question
                shell_block = zeros(size(conv_wm, 1), max_kernel+2);
                shell_block(1:size(conv_wm,1), 1:size(conv_wm,2)) = conv_wm;
                shell_block(1, max_kernel+1) = conv_gm;
                shell_block(1, max_kernel+2) = conv_cs;
                
                % combine the convolution blocks across the shells
                fconv = [fconv; shell_block];
                
            end
            
            % now we do the deconvolution
            % fit SH to each shell
            dwi_sh = [];
            for x=1:nshells
                sh = SH.eval(lmax(x), g_mats{x});
                dwi_sh = [dwi_sh; sh\shells{x}];
            end
            
            % have an initial stab at the fODF fit
            kernel_chop = fconv(1:end, 1:max_kernel);
            f_sh_start = kernel_chop \ dwi_sh;
            n = SH.lmax2n(4);
            f_sh_start(n+1:end,:) = 0; % truncate the higher order components
            f_sh_start = [f_sh_start; 0.1 * ones(2, size(DWI, 2))];
            % few deconvolution parameters
            load dir300.txt
            hr_sh = SH.eval(max(lmax), dir300);
            threshold = par.tau * mean(hr_sh * f_sh_start(1:max_kernel,:), 1);
            lambda = par.lambda * size(fconv, 1) * fconv(1) / size(hr_sh, 1);
            nitts = par.nitts;
            
            % add penalty functions for the isotropic gm and csf components
            line1 = zeros(1, size(hr_sh, 2)+2); line1(end-1) = lambda * hr_sh(1, 1);
            line2 = zeros(1, size(hr_sh, 2)+2); line2( end ) = lambda * hr_sh(1, 1);
            
            hr_sh_combined = [[hr_sh zeros(size(hr_sh, 1), 2)]; line1; line2];
            
            CSD_FOD = zeros(max_kernel+2, size(DWI, 2));
            
            disp('MS deconv...')
            
            parfor x=1:size(dwi_sh, 2)
                CSD_FOD(:, x) = EDTI_Library.E_DTI_MS_deconv_MuSh( f_sh_start(:, x), dwi_sh(:, x),...
                    fconv, hr_sh, hr_sh_combined, threshold(x), lambda, nitts);
            end
            
            disp('MS deconv done!')
            
        end
        
        % From ExploreDTI: helper function for multi-shell deconvolution
        function [ r_sh ] = E_DTI_MS_b0_response_MuSh(dwi, NrB0)
            
            g = zeros(NrB0, 3);
            sh = SH.eval(0, g);
            
            r_sh = zeros(SH.lmax2n(0), 1);
            for i=1:size(dwi, 2)
                r_sh = r_sh + sh\dwi(:,i);
            end
            r_sh = r_sh ./ size(dwi, 2);
        end
        
        % From ExploreDTI: helper function for multi-shell deconvolution
        function [ conv ] = E_DTI_calc_convmat_MuSh( r_sh, target_lmax )
            
            r_sh_n = size(r_sh,1);
            lmax = target_lmax;
            n = SH.lmax2n(lmax);
            if (n > r_sh_n); disp('Error: SH order of response function not high enough for target lmax'); error(''); end;
            d_sh = SH.eval(lmax,[0 0])'; % evaluate at a single Z aligned unit vector ([0 0] = spherical coordinates)
            k = find(d_sh);
            r_rh = r_sh(k)./d_sh(k);
            m = [];
            for l = 0:2:lmax
                m = [m; ones(2*l+1,1)*r_rh(l/2+1)];
            end
            conv = diag(m);
            
            
        end
        
        % From ExploreDTI: helper function for multi-shell deconvolution
        function [f_sh] = E_DTI_MS_deconv_MuSh( f_sh, dwi_sh, fconv, hr_sh, hr_sh_c, threshold, lambda, nitts)
            
            for x=1:nitts
                
                f_hr = hr_sh * f_sh(1:end-2);
                neg  = [f_hr < threshold; f_sh(end-1)<0; f_sh(end)<0];
                
                m = [fconv; lambda*hr_sh_c(neg, :)];
                s = [dwi_sh; zeros(sum(neg), 1)];
                
                f_sh_new = m\s;
                
                if all(f_sh_new == f_sh)
                    f_sh = f_sh_new;
                    break;
                else
                    f_sh = f_sh_new;
                end
            end
            
        end
        
        % From ExploreDTI: helper function to perform MSCSD
        function [r_sh,sel_voxels] = E_DTI_HARDI_CSD_FOD_RC_MuSh_mask(dwi,grad4,bvals,lmax)
            
            % modified from E_DTI_Get_HARDI_CSD_FOD_RC
            % being called by MS_CSD_deconvolution
            % tic
            % load(fin, 'DWI', 'g', 'FA', 'NrB0','b')
            % load(fin,'DWI','VDims','NrB0','b','g','FA','bval','MDims')
            
            % Chantal's recursive calibration code
            % disp('CSD with recursive RF calibration...')
            it = 10;  % maximum iteration times
            
            % try
            %     pctRunOnAll warning off all
            % catch
            %     warning off all
            % end
            
            % h_f = findobj('Tag','MainExploreDTI');
            % data = get(h_f, 'userdata');
            
            t = 0.01; % data.HARDI_GL.RF.RC.PR;
            % lmax = data.HARDI_GL.RF.RC.lmax;
            % suf = data.HARDI_GL.RF.RC.suf;
            
            % DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            %
            % B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, ~isnan(FA), NrB0);
            % DWI(1:NrB0)=B0s;
            % clear B0s;
            %
            % % mask = ~isnan(FA);
            % mask = FA/sqrt(3)>0.01;
            %
            % bvals=round(sum(b(:,[1 4 6]),2)/100)*100;
            % grad=[[zeros(NrB0,3);g],bvals];
            %
            % dwi=vec(E_DTI_DWI_cell2mat(DWI),mask); clear DWI;
            % dwi = single(dwi);
            %
            % M=FA/sqrt(3);
            % M=vec(M,mask);
            %
            % % l=true(size(b,1),1);
            % % dwi = double(dwi(l,:));
            % % grad = grad(l,:);
            %
            % % Select voxels
            % sf_mask=(double(M)>0.01); clear M;
            % % sf_mask = ~isnan(M); clear M;
            %
            % fim = find(sf_mask==1);
            %
            % % if suf>length(fim)/500
            %     suf = round(length(fim)/500);
            % % end
            %
            % fim = fim(1:suf:end);
            % sf_mask(:)=0;
            % sf_mask(fim)=1;
            %
            % EDTI_Library.dwi = EDTI_Library.dwi(:,sf_mask);
            %
            % % select single shell without b0
            % EDTI_Library.dwi = EDTI_Library.dwi(NrB0+1:end,:);
            % grad = grad(NrB0+1:end,1:3);
            
            % dwi in SH
            %%%%%%%%%%%%%%%%%%%% sh = SH(grad,lmax);
            grad = grad4(:,1:3);
            
            sh = SH(lmax,grad);
            
            dwi_sh = sh.coef(dwi);
            
            % axial symmetric RF in z direction for every voxel, FA tensor = 0.05,(0.15)
            % werkt voor lmax = 6 en 8
            
            r_sh = zeros(SH.lmax2n(lmax),size(dwi,2)); % initialization
            
            shcoef = EDTI_Library.E_DTI_create_initial_RF_RC(lmax,max(bvals),0.15,2.1*10^(-3));
            
            for i=1:length(shcoef(:))
                
                r_sh(i,:) = shcoef(i,1);
                
            end
            
            % r_sh(1,:)=1;r_sh(4,:)=-0.0792131438242680;r_sh(11,:)=0.00274230814218417;r_sh(22,:)=-6.26530606143104e-05;
            
            SHPrecomp.init(lmax); % always perform this initialization first
            
            I = cell(1);
            
            for i = 1:it
                %     disp(num2str(i))
                % total RF
                I{i}.nvox = size(r_sh,2); % amount of voxels
                r_shtot = mean(r_sh,2);
                I{i}.r_shtot = r_shtot;
                % amp = sh.amp(r_shtot);
                % plot_odf_Chantal([amp;amp],[grad;-grad]);
                
                % create CSD object
                csd = CSD(r_shtot,lmax); clear r_shtot
                
                % perform constrained deconvolution on SH coefficients
                csd_fod_sh = csd.deconv(dwi_sh);
                
                % peak extraction
                %%%%%%%%%%     [dirs, vals] = SHPrecomp.all_peaks(csd_fod_sh, init_dir, 0, 2);
                [dirs_, vals_] = SHPrecomp.all_peaks(csd_fod_sh, 0, 2);
                
                vals = zeros(2,length(vals_));
                
                for k = 1:length(vals_)
                    if isempty(vals_{k})
                        vals_{k}(1) = nan;
                    end
                    vals(1,k) = vals_{k}(1);
                    if length(vals_{k})>1
                        vals(2,k) = vals_{k}(2);
                    else
                        vals(2,k) = nan;
                    end
                end
                
                dirs = zeros(3,length(vals_));
                
                for k = 1:length(vals_)
                    if isempty(dirs_{k})
                        dirs_{k} = [nan nan nan]';
                    end
                    dirs(1,k) = dirs_{k}(1,1);
                    dirs(2,k) = dirs_{k}(2,1);
                    dirs(3,k) = dirs_{k}(3,1);
                end
                
                I{i}.vals = vals;
                peak_mask = ((vals(2,:)./vals(1,:)<t | isnan(vals(2,:))) & ~isnan(vals(1,:)));
                I{i}.peak_mask = peak_mask;
                dwi_sh = dwi_sh(:,peak_mask);
                r_sh = zeros(SH.lmax2n(lmax),size(dwi_sh,2));
                dirs = dirs(:,peak_mask);
                dwi = dwi(:,peak_mask);
                
                for j = 1:size(dwi_sh,2)
                    fe = dirs(1:3,j);
                    nullsp = null(fe');
                    se = nullsp(:,1);
                    te = nullsp(:,2);
                    rot_grad = zeros(size(grad,1),3);
                    rot = [te se fe];
                    for k = 1:size(grad,1)
                        rot_grad(k,:) = grad(k,:)*rot;
                    end
                    r_sh(:,j) = (SH.eval(lmax,c2s(rot_grad))\dwi(:,j));
                end
                
                j = 0;
                for l_ = 0:2:lmax
                    for m = -l_:l_
                        j = j + 1;
                        if m ~= 0
                            r_sh(j,:) = 0;
                        end
                    end
                end
                
                %     if i==1
                %     save([fin(1:end-4),'_iterations.mat'],'I');
                %     disp('saving iterations...')
                %     end
                %     if i>1
                %         if (abs(I{i}.r_shtot-I{i-1}.r_shtot)./I{i-1}.r_shtot)<0.01
                %             break
                %         end
                %     end
                
                if i>1
                    %        save([fin(1:end-4),'_iterations.mat'],'I','-append');
                    change = abs(I{i}.r_shtot-I{i-1}.r_shtot)./I{i-1}.r_shtot;
                    if  all(change(~isnan(change))<0.01)
                        break
                    end
                end
                
                
            end
            %%%%%%%% save I %%%%
            
            % build sf_mask
            sel_voxels = I{i}.peak_mask;
            for It = length(I)-1:-1:1
                new_mask = zeros(I{It}.nvox,1);
                new_mask(I{It}.peak_mask) = sel_voxels;
                sel_voxels = new_mask;
            end
            
            % save([p1 filesep p2 '_I.mat'],'I');
            % disp('saving I...')
            
            % create sf_mask
            
            % masks = zeros(size(I,2),I{1,1}.nvox);
            %
            % five_masks = zeros(size(I,2),length(find(sf_mask)));
            %
            %
            % for i = 1:size(I,2)
            %     m = ones(1,I{i}.nvox);
            %     for j = (i-1):-1:1
            %         m = unvec_Chantal(m,I{j}.peak_mask);
            %     end
            %     m = ~isnan(m);
            %     masks(i,m) = 1;
            % end
            %
            % % %% sf_mask = sum(unvec(masks,mask),4); %save([directory,tag,'mask.mat'],'ma');
            % five_masks(:,1:suf:end) = masks(:,:);
            %
            % sf = unvec(five_masks,mask);
            % % sf = unvec(masks,mask);
            % sf_mask_end = sf(:,:,:,length(I));
            
            % % sf_mask(isnan(sf_mask)) = 0;
            % save([p1 filesep p2 '_sf_mask.mat'],'sf_mask_end','sf');
            % disp('saving sf_mask...')
            
            % EDTI_Library.E_DTI_write_nifti_file(sf,VDims,fn_sf);
            
            % Plot mean PR and nr of voxels over iteration
            % % nvox=[];
            % % ratiovals=[];
            
            r_shtot=[];
            
            for i = 1:size(I,2)
                %     nvox=[nvox,I{i}.nvox];
                %     ratiovals=[ratiovals,mean([I{i}.vals(2,~isnan(I{i}.vals(2,:)))/I{i}.vals(1,~isnan(I{i}.vals(2,:))),zeros(1,length(I{i}.vals(2,isnan(I{i}.vals(2,:)))))])];
                r_shtot = [r_shtot,I{i}.r_shtot];
            end
            
            % figure
            % plot(nvox), xlabel Iteration, ylabel('Number of voxels used for RF calculation'); saveas(gcf,[directory,tag,'_nvoxels.png']);
            % figure
            % plot(ratiovals), xlabel Iteration, ylabel('Mean ratio magnitude second peak / first peak'); saveas(gcf,[directory,tag,'_meanratio.png']);
            %
            % f = load(['icosahedron4.mat']);
            % g2 = f.Expression1;
            % %%% sh2 = SH(g2, lmax);
            % sh2 = SH(lmax, g2);
            
            % Load DTI.mat file
            % load([directory,file],'DWI','VDims','b','g','NrB0','MDims','FEFA','FA','FE','SE','eigval','DT');
            % load(fin,'DWI','VDims','NrB0','b','g','FA','bval','MDims')
            
            % % bvals = sum(b(:,[1 4 6]),2);
            % bvals = round(sum(b(:,[1 4 6]),2)/100)*100;
            % grad = [[zeros(NrB0,3);g],bvals];
            % dwi = vec(E_DTI_DWI_cell2mat(DWI),mask); clear DWI;
            % dwi = single(dwi);
            
            % Do CSD for iteration(s) iter, mostly last iteration
            iter = size(I,2);
            
            % % select single shell
            % EDTI_Library.dwi = EDTI_Library.dwi(NrB0+1:end,:);
            % % b = b(NrB0+1:end,:);
            % grad = grad(NrB0+1:end,1:3);
            
            % create SH object
            % % sh = SH(grad, lmax);
            % sh = SH(lmax, grad);
            %
            % % estimate SH coefficients
            % dwi_sh = sh.coef(dwi);
            
            for i = iter
                %     disp(num2str(i))
                r_sh = r_shtot(:,i);
            end
            
            %%%%%%% save RF to a folder
            
            % if ~exist([p1 '_RF'], 'dir')
            %   mkdir([p1 '_RF']);
            % end
            %
            % dirname_RF = [p1 '_RF'];
            %
            % save([dirname_RF filesep p2 '_RF.mat'],'r_sh');
            % disp('saving RF...')
            % toc
            
        end
        
        % Helper function for multi-shell deconvolution
        function [mask,s_mask] = UniformMaskSubsamp(mask,desired_points)
            wm_sampling_ratio = size(mask,2)/desired_points;
            if(wm_sampling_ratio < 1)
                wm_sampling_ratio = 1;
            else
                wm_sampling_ratio = round(wm_sampling_ratio);
            end
            s_mask = 1:size(mask,2);
            s_mask = s_mask(1:wm_sampling_ratio:end);
            mask = mask(:, 1:wm_sampling_ratio:end);
        end
        
        % From ExploreDTI: Perform fiber tractography of any FOD in SH - Adapted
        function WholeBrainFODTractography(reference_mat,CSD_FOD_or_file,p,f_out)
            global MRIToolkit
            
            tic
            
            f_in = reference_mat;
            
            load(f_in,'VDims')
            
            if(ischar(CSD_FOD_or_file))
                CSD_FOD = EDTI_Library.E_DTI_read_nifti_file(CSD_FOD_or_file);
            else
                CSD_FOD = CSD_FOD_or_file;
                clear CSD_FOD_or_file;
            end
            
            SR = p.SeedPointRes;
            mask_s = ~isnan(CSD_FOD(:,:,:,1));
            
            disp('Calculating seed points...')
            if(isfield(p,'SeedMask') && ~isempty(p.SeedMask))
                if(ischar(p.SeedMask))
                    seed_mask = EDTI_Library.E_DTI_read_nifti_file(p.SeedMask);
                    p.SeedMask = seed_mask;
                end
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(p.SeedMask > 0, SR, VDims, p);
            else
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(mask_s, SR, VDims, p);
            end
            
            if isempty(seedPoint)
                disp('No seed points found - processing stopped.')
            else
                disp('Seed point calculation done.')
            end
            
            disp('Calculating trajectories...')
            
            stepSize = p.StepSize;
            threshold = p.blob_T;
            maxAngle = p.AngleThresh;
            lengthRange = [p.FiberLengthRange(1) p.FiberLengthRange(2)];
            v2w = diag(VDims); v2w(4,4) = 1;
            
            t = [];
            if(isfield(MRIToolkit,'fibertracker') < 1 || isfield(MRIToolkit.fibertracker,'type') < 1 ...
                    || isempty(MRIToolkit.fibertracker.type))
                t = SHTracker(v2w);
            else
                eval(['t = ' MRIToolkit.fibertracker.type '(v2w);']);
                if(isfield(MRIToolkit.fibertracker,'parameters'))
                    fields = fieldnames(MRIToolkit.fibertracker.parameters);
                    for field_id=1:length(fields)
                        eval(['t.' fields{field_id} '= MRIToolkit.fibertracker.parameters.' fields{field_id} ';']);
                    end
                end
            end
            t.setData(CSD_FOD);
            t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(false);
            [Tracts, TractsFOD] = t.track(seedPoint);
            
            num_tracts = size(Tracts,2);
            TractL = cell(1,num_tracts);
            for i = 1:num_tracts
                TractL{i} = single(size(Tracts{i},1));
            end
            
            TL = cell2mat(TractL(:))*stepSize;
            IND = and(TL>=lengthRange(1),TL<=lengthRange(2));
            Tracts = Tracts(IND);
            TractsFOD = TractsFOD(IND);
            
            num_tracts = size(Tracts,2);
            TractL = cell(1,num_tracts);
            for i = 1:num_tracts
                TractL{i} = single(size(Tracts{i},1));
            end
            
            disp('Trajectory calculations done.')
            
            disp('Processing diffusion info along the tracts...');
            
            load(f_in,'DT','MDims');
            
            [TractFA, TractFE, TractAng, TractGEO, TractLambdas, TractMD] =...
                EDTI_Library.E_DTI_diff_measures_vectorized(Tracts, VDims, DT);
            
            % TractL = cellfun(@length, Tracts, 'UniformOutput', 0);
            
            FList = (1:length(Tracts))';
            
            TractMask = repmat(0,MDims);
            
            for i = 1:length(FList)
                IND = unique(sub2ind(MDims,...
                    round(double(Tracts{i}(:,1))./VDims(1)),...
                    round(double(Tracts{i}(:,2))./VDims(2)),...
                    round(double(Tracts{i}(:,3))./VDims(3))));
                TractMask(IND) = TractMask(IND) + 1;
            end
            
            disp('Processing diffusion info along the tracts done.');
            
            Parameters.Step_size = stepSize;
            Parameters.FOD_threshold = threshold;
            Parameters.Angle_threshold = maxAngle;
            Parameters.Length_range = lengthRange;
            Parameters.Seed_resolution = SR;
            
            
            disp('Saving trajectories...')
            save(f_out,'FList','Tracts','TractL','TractFE',...
                'TractFA','TractAng','TractGEO','TractMD',...
                'TractLambdas','TractMask','VDims','TractsFOD','Parameters');
            disp('Saving trajectories done.')
            
            ti = toc;
            
            m=ti/60;
            disp(['Tracking (CSD - FOD interpolation) computation time was ' num2str(m) ' min.'])
            
        end
        
        % From ExploreDTI: List of possible exports
        function [List_var, LE] = E_DTI_Complete_List_var
            
            % HE = E_DTI_Get_HARDI_ext;
            HE = '_HARDI.nii';
            
            c=0;
            c=c+1;List_var{c} = 'Fractional anisotropy (''_FA.nii'')'; %1
            c=c+1;List_var{c} = 'Mean diffusivity (''_MD.nii'')'; %2
            c=c+1;List_var{c} = 'Largest eigenvalue (''_L1.nii'')'; %3
            c=c+1;List_var{c} = 'Middle eigenvalue (''_L2.nii'')'; %4
            c=c+1;List_var{c} = 'Smallest eigenvalue (''_L3.nii'')'; %5
            c=c+1;List_var{c} = 'Radial/transverse diffusivity (''_RD.nii'')'; %6
            c=c+1;List_var{c} = 'Relative Anisotropy (''_RA.nii'')'; %7
            c=c+1;List_var{c} = 'Estimated B0 image (''_B0.nii'')'; %8
            c=c+1;List_var{c} = 'Diffusion tensor (''_DT.nii'')'; %9
            c=c+1;List_var{c} = 'Westin-measure CL (''_CL.nii'')'; %10
            c=c+1;List_var{c} = 'Westin-measure CP (''_CP.nii'')'; %11
            c=c+1;List_var{c} = 'Westin-measure CS (''_CS.nii'')'; %12
            c=c+1;List_var{c} = 'First eigenvector (''_FE.nii'')'; %13
            c=c+1;List_var{c} = 'Second eigenvector (''_SE.nii'')'; %14
            c=c+1;List_var{c} = 'Third eigenvector (''_TE.nii'')'; %15
            c=c+1;List_var{c} = '|first eigenvector| (''_abs_FE.nii'')'; %16
            c=c+1;List_var{c} = '|second eigenvector| (''_abs_SE.nii'')'; %17
            c=c+1;List_var{c} = '|third eigenvector| (''_abs_TE.nii'')'; %18
            c=c+1;List_var{c} = 'FA x |first eigenvector| (''_FA_abs_FE.nii'')'; %19
            c=c+1;List_var{c} = 'CL x |first eigenvector| (''_CL_abs_FE.nii'')'; %20
            c=c+1;List_var{c} = 'CP x |third eigenvector| (''_CP_abs_TE.nii'')'; %21
            c=c+1;List_var{c} = 'Mean of DWI residuals from DT (''_mean_res_DWI_DT.nii'')'; %22
            c=c+1;List_var{c} = 'Max of DWI residuals from DT (''_max_res_DWI_DT.nii'')'; %23
            c=c+1;List_var{c} = 'Mean of B0 residuals from DT (''_mean_res_B0_DT.nii'')'; %24
            c=c+1;List_var{c} = 'Max of B0 residuals from DT (''_max_res_B0_DT.nii'')'; %25
            c=c+1;List_var{c} = 'Standard deviation DWIs (''_std_DWI.nii'')'; %26
            c=c+1;List_var{c} = 'Standard deviation B0s (''_std_B0.nii'')'; %27
            c=c+1;List_var{c} = 'Skew eigenvalues (''_skew_L.nii'')'; %28
            c=c+1;List_var{c} = 'DWI>B0 values on FA (''_PNPIV.nii'')'; %29
            c=c+1;List_var{c} = 'DWIs with B0(s) (''_DWIs.nii'')'; %30
            c=c+1;List_var{c} = 'Local dyadic coherence (''_Kappa.nii'')'; %31
            c=c+1;List_var{c} = 'Inter-voxel diffusion coherence (''_IVDC.nii'')'; %32
            c=c+1;List_var{c} = 'B-matrix (''_B_matrix.txt'')'; %33
            c=c+1;List_var{c} = 'Gradient directions (''_grads.txt'')'; %34
            c=c+1;List_var{c} = 'Apparent diffusion coefficients (''_ADCs.nii'')'; %35
            c=c+1;List_var{c} = 'ADCs - reconstructed with DT model (''_ADCs_DT.nii'')'; %36
            c=c+1;List_var{c} = 'DWIs - reconstructed with DT model (''_DWIs_DT.nii'')'; %37
            c=c+1;List_var{c} = 'DWI residuals from DT (''_res_DWI_DT.nii'')'; %38
            c=c+1;List_var{c} = 'ADC residuals from DT (''_res_ADC_DT.nii'')'; %39
            c=c+1;List_var{c} = 'Mask (''_mask.nii'')'; %40
            c=c+1;List_var{c} = 'Lattice Index (''_LI.nii'')'; %41
            c=c+1;List_var{c} = 'Mean of DWIs (''_mean_DWIs.nii'')'; %42
            c=c+1;List_var{c} = 'Mean of B0s (''_mean_B0s.nii'')'; %43
            c=c+1;List_var{c} = ['HARDI (see settings) (''' HE ''')']; %44
            c=c+1;List_var{c} = 'Diffusion kurtosis tensor (''_KT.nii'')'; %45
            c=c+1;List_var{c} = 'Mean kurtosis (''_MK.nii'')'; %46
            c=c+1;List_var{c} = 'Axial kurtosis (''_AK.nii'')'; %47
            c=c+1;List_var{c} = 'Radial kurtosis (''_RK.nii'')'; %48
            c=c+1;List_var{c} = 'Kurtosis anisotropy (''_KA.nii'')'; %49
            c=c+1;List_var{c} = 'Mode of anisotropy (''_AM.nii'')'; %50
            c=c+1;List_var{c} = 'Diffusion tensor norm (''_DTN.nii'')'; %51
            c=c+1;List_var{c} = 'Skeletonized FA (''_FA_Skeleton.nii'')'; %52
            c=c+1;List_var{c} = 'Axonal Water Fraction (from DKI) (''_AWF.nii'')'; %53
            c=c+1;List_var{c} = 'Tortuosity (from DKI) (''_TORT.nii'')'; %54
            c=c+1;List_var{c} = 'Axial extra-axonal diffusivity (from DKI) (''_AxEAD.nii'')'; %55
            c=c+1;List_var{c} = 'Radial extra-axonal diffusivity (from DKI) (''_RadEAD.nii'')'; %56
            
            c=0;
            c=c+1;LE{c} = '_FA.nii';
            c=c+1;LE{c} = '_MD.nii';
            c=c+1;LE{c} = '_L1.nii';
            c=c+1;LE{c} = '_L2.nii';
            c=c+1;LE{c} = '_L3.nii';
            c=c+1;LE{c} = '_RD.nii';
            c=c+1;LE{c} = '_RA.nii';
            c=c+1;LE{c} = '_B0.nii';
            c=c+1;LE{c} = '_DT.nii';
            c=c+1;LE{c} = '_CL.nii';
            c=c+1;LE{c} = '_CP.nii';
            c=c+1;LE{c} = '_CS.nii';
            c=c+1;LE{c} = '_FE.nii';
            c=c+1;LE{c} = '_SE.nii';
            c=c+1;LE{c} = '_TE.nii';
            c=c+1;LE{c} = '_abs_FE.nii';
            c=c+1;LE{c} = '_abs_SE.nii';
            c=c+1;LE{c} = '_abs_TE.nii';
            c=c+1;LE{c} = '_FA_abs_FE.nii';
            c=c+1;LE{c} = '_CL_abs_FE.nii';
            c=c+1;LE{c} = '_CP_abs_TE.nii';
            c=c+1;LE{c} = '_mean_res_DWI_DT.nii';
            c=c+1;LE{c} = '_max_res_DWI_DT.nii';
            c=c+1;LE{c} = '_mean_res_B0_DT.nii';
            c=c+1;LE{c} = '_max_res_B0_DT.nii';
            c=c+1;LE{c} = '_std_DWI.nii';
            c=c+1;LE{c} = '_std_B0.nii';
            c=c+1;LE{c} = '_skew_L.nii';
            c=c+1;LE{c} = '_PNPIV.nii';
            c=c+1;LE{c} = '_DWIs.nii';
            c=c+1;LE{c} = '_Kappa.nii';
            c=c+1;LE{c} = '_IVDC.nii';
            c=c+1;LE{c} = '_B_matrix.txt';
            c=c+1;LE{c} = '_grads.txt';
            c=c+1;LE{c} = '_ADCs.nii';
            c=c+1;LE{c} = '_ADCs_DT.nii';
            c=c+1;LE{c} = '_DWIs_DT.nii';
            c=c+1;LE{c} = '_res_DWI_DT.nii';
            c=c+1;LE{c} = '_res_ADC_DT.nii';
            c=c+1;LE{c} = '_mask.nii';
            c=c+1;LE{c} = '_LI.nii';
            c=c+1;LE{c} = '_mean_DWIs.nii';
            c=c+1;LE{c} = '_mean_B0s.nii';
            c=c+1;LE{c} = '_hardi.nii';%HE;
            c=c+1;LE{c} = '_KT.nii';
            c=c+1;LE{c} = '_MK.nii';
            c=c+1;LE{c} = '_AK.nii';
            c=c+1;LE{c} = '_RK.nii';
            c=c+1;LE{c} = '_KA.nii';
            c=c+1;LE{c} = '_AM.nii';
            c=c+1;LE{c} = '_DTN.nii';
            c=c+1;LE{c} = '_FA_Skeleton.nii';
            c=c+1;LE{c} = '_AWF.nii';
            c=c+1;LE{c} = '_TORT.nii';
            c=c+1;LE{c} = '_AxEAD.nii';
            c=c+1;LE{c} = '_RadEAD.nii';
            
        end
        
        % From ExploreDTI: helper for kurtosis export
        function fa = E_DTI_FrAn(L1,L2,L3)
            
            mask = ~isnan(L1);
            L1 = L1(mask);
            L2 = L2(mask);
            L3 = L3(mask);
            
            D1 = L1-L2;
            D1 = D1.^2;
            
            D2 = L2-L3;
            D2 = D2.^2;
            
            D1 = D1 + D2;
            
            D2 = L3-L1;
            D2 = D2.^2;
            
            D1 = D1 + D2;
            D1 = sqrt(D1);
            
            clear D2;
            
            L1 = L1.^2;
            L2 = L2.^2;
            L3 = L3.^2;
            
            L1 = L1 + L2;
            L1 = L1 + L3;
            
            clear L2 L3;
            
            L1 = 2*L1;
            L1 = sqrt(L1);
            fa = repmat(double(nan),size(mask));
            
            warning off all
            dummy = D1./L1;
            dummy(isnan(dummy))=0;
            dummy(isinf(dummy))=0;
            fa(mask) = dummy;
            warning on all
        end
        
        % From ExploreDTI: helper for kurtosis export
        function [FEFA, FA, FE, SE, eigval] = E_DTI_View3DDTMaps_vec(DT, mask)
            
            if iscell(DT)
                MDims = size(DT{1});
            else
                MDims = size(DT(:,:,:,1));
            end
            
            Dummy = repmat(double(nan),MDims);
            
            try
                X = repmat(double(0),[6 sum(mask(:))]);
            catch
                X = repmat(single(0),[6 sum(mask(:))]);
            end
            
            if iscell(DT)
                
                for k=1:size(X,1);
                    Dummy = DT{k};
                    X(k,:)=double(Dummy(mask(:)));
                end
                
                
            else
                
                order = [1 2 3 5 6 9];
                
                for k=1:size(X,1);
                    Dummy = DT(:,:,:,order(k));
                    X(k,:)=double(Dummy(mask(:)));
                end
                
            end
            
            % [V,E]=eig_Lex(double(X));
            
            [V, E, PD] = EDTI_Library.E_DTI_eig_Lex_H(X);
            clear X;
            
            % PosDef = repmat(logical(0), [size(DT,1) size(DT,2) size(DT,3)]);
            % PosDef(mask)=~PD;
            % imagescn(PosDef,[0 1])
            try
                FE = repmat(double(nan), [MDims 3]);
                SE = repmat(double(nan), [MDims 3]);
                eigval = repmat(double(nan), [MDims 3]);
            catch
                FE = repmat(single(nan), [MDims 3]);
                SE = repmat(single(nan), [MDims 3]);
                eigval = repmat(single(nan), [MDims 3]);
            end
            
            for i=1:3
                Dummy(mask) = E(i,:);
                eigval(:,:,:,i) = Dummy;
                Dummy(mask) = V(i,1,:);
                FE(:,:,:,i) = Dummy;
                Dummy(mask) = V(i,2,:);
                SE(:,:,:,i) = Dummy;
            end
            clear V;
            
            try
                FA = repmat(double(nan), MDims);
            catch
                FA = repmat(single(nan), MDims);
            end
            
            FA(mask) = EDTI_Library.E_DTI_FrAn(E(1,:),E(2,:),E(3,:));
            
            clear E;
            try
                FEFA = repmat(double(nan), [MDims 3]);
            catch
                FEFA = repmat(single(nan), [MDims 3]);
            end
            
            for i=1:3
                Dummy = FE(:,:,:,i);
                Dummy(mask) = FA(mask).*abs(Dummy(mask));
                FEFA(:,:,:,i) = Dummy;
            end
            FEFA(FEFA<0)=0;
            FEFA(FEFA>1)=1;
            FA(mask) = FA(mask)*sqrt(3);
        end
        
        % From ExploreDTI: helper for kurtosis export
        function AK = E_DTI_Axial_Kurtosis_c(KT,DT)
            
            DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            mask = ~isnan(KT{1});
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_View3DDTMaps_vec(DT, mask);
            
            f1 = FE(:,:,:,1);
            f2 = FE(:,:,:,2);
            f3 = FE(:,:,:,3);
            
            f1 = f1(mask);
            f2 = f2(mask);
            f3 = f3(mask);
            
            g = [f1 f2 f3];
            
            clear FEFA FA FE SE eigval;
            
            MDsq = repmat(single(nan),size(KT{1}));
            AK = repmat(single(nan),size(KT{1}));
            AK(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i});
            end
            
            dt = [DT{1} DT{4} DT{6} DT{2} DT{3} DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i});
            end
            
            for i = [1 2 3 10 11 12]
                %     KT{i}(KT{i}<=0) = 0;
                KT{i}(KT{i}<=0) = abs(KT{i}(KT{i}<=0));
            end
            
            kt = [KT{1} KT{2} KT{3} KT{4} KT{5} KT{6} ...
                KT{7} KT{8} KT{9} KT{10} KT{11} KT{12} ...
                KT{13} KT{14} KT{15}];
            
            clear KT;
            
            
            dum_adc = (sum(B.*dt,2)).^2;
            dum_kt = sum(A.*kt,2);
            
            AK(mask) = (MDsq(mask).*dum_kt(:))./dum_adc(:);
            
            [dummy, AK] = EDTI_Library.E_DTI_local_outlier_correction_avg(AK,mask);
            % y = prctile(AK(mask),[1 99]);
            % AK = EDTI_Library.E_DTI_outlier_smoothing(AK,mask,y);
        end
        
        % From ExploreDTI: helper for kurtosis export
        function AK = E_DTI_Axial_Kurtosis(KT,DT)
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            global MRIToolkit;
            if(isfield(MRIToolkit,'DKI_cleanup') && MRIToolkit.DKI_cleanup == true)
                AK = EDTI_Library.E_DTI_Axial_Kurtosis_c(KT,DT);
                return;
            end
            
            mask = ~isnan(KT{1});
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_View3DDTMaps_vec(DT, mask);
            
            f1 = FE(:,:,:,1);
            f2 = FE(:,:,:,2);
            f3 = FE(:,:,:,3);
            
            f1 = f1(mask);
            f2 = f2(mask);
            f3 = f3(mask);
            
            g = [f1 f2 f3];
            
            clear FEFA FA FE SE eigval;
            
            MDsq = repmat(single(nan),size(KT{1}));
            AK = repmat(single(nan),size(KT{1}));
            AK(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i});
            end
            
            dt = [DT{1} DT{4} DT{6} DT{2} DT{3} DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i});
            end
            
            kt = [KT{1} KT{2} KT{3} KT{4} KT{5} KT{6} ...
                KT{7} KT{8} KT{9} KT{10} KT{11} KT{12} ...
                KT{13} KT{14} KT{15}];
            
            clear KT;
            
            % for i=1:size(g,1)
            
            dum_adc = (sum(B.*dt,2)).^2;
            dum_kt = sum(A.*kt,2);
            
            AK(mask) = (MDsq(mask).*dum_kt(:))./dum_adc(:);
            
            % y = prctile(AK(mask),[1 99]);
            %
            % AK = EDTI_Library.E_DTI_outlier_smoothing(AK,mask,y);
        end
        
        % From ExploreDTI: helper for kurtosis export
        function RK = E_DTI_Radial_Kurtosis_c(KT,DT)
            
            DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            N=10;
            
            mask = ~isnan(KT{1});
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_View3DDTMaps_vec(DT, mask);
            
            TE = cross(FE,SE,4);
            
            s1 = SE(:,:,:,1);
            s2 = SE(:,:,:,2);
            s3 = SE(:,:,:,3);
            
            t1 = TE(:,:,:,1);
            t2 = TE(:,:,:,2);
            t3 = TE(:,:,:,3);
            
            s1 = s1(mask);
            s2 = s2(mask);
            s3 = s3(mask);
            
            t1 = t1(mask);
            t2 = t2(mask);
            t3 = t3(mask);
            
            clear FEFA FA FE SE TE eigval;
            
            
            MDsq = repmat(single(nan),size(KT{1}));
            RK = repmat(single(nan),size(KT{1}));
            RK(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i});
            end
            
            dt = [DT{1} DT{4} DT{6} DT{2} DT{3} DT{5}];
            clear DT;
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i});
            end
            
            for i = [1 2 3 10 11 12]
                %     KT{i}(KT{i}<=0) = 0;
                KT{i}(KT{i}<=0) = abs(KT{i}(KT{i}<=0));
            end
            
            kt = [KT{1} KT{2} KT{3} KT{4} KT{5} KT{6} ...
                KT{7} KT{8} KT{9} KT{10} KT{11} KT{12} ...
                KT{13} KT{14} KT{15}];
            clear KT;
            
            th = linspace(0,2*pi,N);
            
            for i=1:N-1
                
                g = [s1 s2 s3]*cos(th(i)) + [t1 t2 t3]*sin(th(i));
                
                
                A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
                
                
                dum_adc = (sum(B.*dt,2)).^2;
                dum_kt = sum(A.*kt,2);
                
                RK(mask) = RK(mask) + (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
            end
            
            RK(mask) = RK(mask)/(N-1);
            
            [dummy, RK] = EDTI_Library.E_DTI_local_outlier_correction_avg(RK,mask);
            % y = prctile(RK(mask),[1 99]);
            % RK = EDTI_Library.E_DTI_outlier_smoothing(RK,mask,y);
        end
        
        % From ExploreDTI: helper for kurtosis export
        function RK = E_DTI_Radial_Kurtosis(KT,DT)
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            global MRIToolkit;
            if(isfield(MRIToolkit,'DKI_cleanup') && MRIToolkit.DKI_cleanup == true)
                RK = EDTI_Library.E_DTI_Radial_Kurtosis_c(KT,DT);
                return;
            end
            
            N = 10;
            
            mask = ~isnan(KT{1});
            
            [FEFA, ~, FE, SE, eigval] = EDTI_Library.E_DTI_View3DDTMaps_vec(DT, mask);
            
            TE = cross(FE,SE,4);
            
            s1 = SE(:,:,:,1);
            s2 = SE(:,:,:,2);
            s3 = SE(:,:,:,3);
            
            t1 = TE(:,:,:,1);
            t2 = TE(:,:,:,2);
            t3 = TE(:,:,:,3);
            
            s1 = s1(mask);
            s2 = s2(mask);
            s3 = s3(mask);
            
            t1 = t1(mask);
            t2 = t2(mask);
            t3 = t3(mask);
            
            clear FEFA FA FE SE TE eigval;
            
            
            MDsq = repmat(single(nan),size(KT{1}));
            RK = repmat(single(nan),size(KT{1}));
            RK(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i});
            end
            
            dt = [DT{1} DT{4} DT{6} DT{2} DT{3} DT{5}];
            clear DT;
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i});
            end
            
            kt = [KT{1} KT{2} KT{3} KT{4} KT{5} KT{6} ...
                KT{7} KT{8} KT{9} KT{10} KT{11} KT{12} ...
                KT{13} KT{14} KT{15}];
            clear KT;
            
            th = linspace(0,2*pi,N);
            
            for i=1:N-1
                
                g = [s1 s2 s3]*cos(th(i)) + [t1 t2 t3]*sin(th(i));
                
                
                A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                    4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                    4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                    4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                    6*(g(:,1).^2).*(g(:,2).^2) ...
                    6*(g(:,1).^2).*(g(:,3).^2) ...
                    6*(g(:,2).^2).*(g(:,3).^2) ...
                    12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                    12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                    12*g(:,1).*g(:,2).*(g(:,3).^2)];
                
                B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
                
                
                dum_adc = (sum(B.*dt,2)).^2;
                dum_kt = sum(A.*kt,2);
                
                RK(mask) = RK(mask) + (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
            end
            
            RK(mask) = RK(mask)/(N-1);
            
            % y = prctile(RK(mask),[1 99]);
            % RK = EDTI_Library.E_DTI_outlier_smoothing(RK,mask,y);
        end
        
        % From ExploreDTI: helper for kurtosis export
        function KA = E_DTI_Kurtosis_Anisotropy_c(KT,DT)
            
            mask = ~isnan(KT{1});
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i})';
            end
            
            for i=[1 4 6]
                %     DT{i}(DT{i}<=0) = 0;
                DT{i}(DT{i}<0) = abs(DT{i}(DT{i}<0));
            end
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            p = which('Grad_dirs_1024.txt');
            g = textread(p);
            
            
            MDsq = repmat(single(nan),size(KT{1}));
            MK = repmat(single(nan),size(KT{1}));
            MK(mask)=0;
            KA = repmat(single(nan),size(KT{1}));
            KA(mask)=0;
            
            MDsq(mask) = ((DT{1}+DT{4}+DT{6})/3).^2;
            
            dt = [DT{1}; DT{4}; DT{6}; DT{2}; DT{3}; DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i})';
            end
            
            %-
            for i = [1 2 3 10 11 12]
                %     KT{i}(KT{i}<=0) = 0;
                KT{i}(KT{i}<=0) = abs(KT{i}(KT{i}<=0));
            end
            
            kt = [KT{1}; KT{2}; KT{3}; KT{4}; KT{5}; KT{6}; ...
                KT{7}; KT{8}; KT{9}; KT{10}; KT{11}; KT{12}; ...
                KT{13}; KT{14}; KT{15}];
            
            clear KT;
            
            cn = MDsq(mask);
            cn(:) = 0;
            V_MK = cn;
            V_KA = cn;
            VMDsq = MDsq(mask);
            
            for i=1:size(g,1)
                
                dum_adc = (B(i,:)*dt).^2;
                dum_kt = A(i,:)*kt;
                
                dum_adc = dum_adc(:);
                dum_kt = dum_kt(:);
                
                M_c = dum_kt(:)>0;
                cn(M_c) = cn(M_c)+1;
                
                %     MK(mask(M_c)) = MK(mask) + (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
                dum = (VMDsq(M_c).*dum_kt(M_c))./dum_adc(M_c);
                
                V_MK(M_c) = V_MK(M_c) + dum;
                V_KA(M_c) = V_KA(M_c) + dum.^2;
            end
            
            V_MK = V_MK./cn;
            V_KA = V_KA./cn;
            
            V_MK = V_MK.^2;
            
            V_KA = sqrt(V_KA-V_MK);
            
            KA(mask) = V_KA(:);
            KA = real(KA);

            [dummy, KA] = EDTI_Library.E_DTI_local_outlier_correction_avg(KA,mask);
            % y = prctile(KA(mask),[0.1 99.9]);
            % KA = EDTI_Library.E_DTI_outlier_smoothing(KA,mask,y);
            
        end
        
        % From ExploreDTI: helper for kurtosis export
        function KA = E_DTI_Kurtosis_Anisotropy(KT,DT)
            
            global MRIToolkit;
            if(isfield(MRIToolkit,'DKI_cleanup') && MRIToolkit.DKI_cleanup == true)
                KA = EDTI_Library.E_DTI_Kurtosis_Anisotropy_c(KT,DT);
                return;
            end
            
            % DT{1} = abs(DT{1});DT{4} = abs(DT{4});DT{6} = abs(DT{6});
            
            p = which('Grad_dirs_1024.txt');
            g = textread(p);
            
            mask = ~isnan(KT{1});
            
            MDsq = repmat(single(nan),size(KT{1}));
            MK = repmat(single(nan),size(KT{1}));
            MK(mask)=0;
            KA = repmat(single(nan),size(KT{1}));
            KA(mask)=0;
            MDsq(mask) = ((DT{1}(mask)+DT{4}(mask)+DT{6}(mask))/3).^2;
            
            for i=1:6
                DT{i} = DT{i}(mask);
                DT{i} = double(DT{i})';
            end
            
            dt = [DT{1}; DT{4}; DT{6}; DT{2}; DT{3}; DT{5}];
            
            clear DT;
            
            A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
                4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
                4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
                4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
                6*(g(:,1).^2).*(g(:,2).^2) ...
                6*(g(:,1).^2).*(g(:,3).^2) ...
                6*(g(:,2).^2).*(g(:,3).^2) ...
                12*g(:,2).*g(:,3).*(g(:,1).^2) ...
                12*g(:,1).*g(:,3).*(g(:,2).^2) ...
                12*g(:,1).*g(:,2).*(g(:,3).^2)];
            
            B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];
            
            
            for i=1:15
                KT{i} = KT{i}(mask);
                KT{i} = double(KT{i})';
            end
            
            kt = [KT{1}; KT{2}; KT{3}; KT{4}; KT{5}; KT{6}; ...
                KT{7}; KT{8}; KT{9}; KT{10}; KT{11}; KT{12}; ...
                KT{13}; KT{14}; KT{15}];
            
            clear KT;
            
            for i=1:size(g,1)
                
                dum_adc = (B(i,:)*dt).^2;
                dum_kt = A(i,:)*kt;
                
                dum = (MDsq(mask).*dum_kt(:))./dum_adc(:);
                
                MK(mask) = MK(mask) + dum;
                KA(mask) = KA(mask) + dum.^2;
                
            end
            
            MK(mask) = MK(mask)/size(g,1);
            KA(mask) = KA(mask)/size(g,1);
            
            MK(mask) = MK(mask).^2;
            
            KA(mask) = sqrt(KA(mask)-MK(mask));
            
            % y = prctile(KA(mask),[1 99]);
            % KA = EDTI_Library.E_DTI_outlier_smoothing(KA,mask,y);
        end
        
        % From ExploreDTI: Export the metrics stored in the MAT file to various
        % formats
        function E_DTI_Convert_mat_2_nii(fin,folout,LS)
            
            disp(['Converting stuff from ''' fin ''' to *.nii files ...'])
            
            [PATHSTR,NAME,EXT] = fileparts(fin);
            
            warning off all
            
            load(fin,'VDims')
            
            if ~exist('VDims','var')
                disp(['Format of ''' fin ''' not correct, skipping data!'])
                return;
            end
            
            [CLL, LE] = EDTI_Library.E_DTI_Complete_List_var;
            
            for i=1:length(LS)
                data_v{i} = [];
                fns{i} = [];
            end
            
            for i=1:length(LS)
                
                if strcmp(LS{i},CLL{1}) %FA
                    
                    %         load(fin,'FA')
                    %
                    %         if ~exist('FA','var')
                    %             disp(['Format of ''' fin ''' not correct, skipping data!'])
                    %             continue;
                    %         end
                    %
                    %         FA(isnan(FA))=0;
                    %         FA = FA/sqrt(3);
                    %         FA(FA>1)=1;
                    %         FA(FA<0)=0;
                    %         data_v{i} = FA;
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = EDTI_Library.FrAn_calc(eigval(:,:,:,1),eigval(:,:,:,2),eigval(:,:,:,3));
                    
                    fns{i} = LE{1};
                    
                elseif strcmp(LS{i},CLL{2}) %MD
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    data_v{i} = mean(eigval,4);
                    fns{i} = LE{2};
                    
                elseif strcmp(LS{i},CLL{3}) %L1
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = eigval(:,:,:,1);
                    fns{i} = LE{3};
                    
                elseif strcmp(LS{i},CLL{4}) %L2
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = eigval(:,:,:,2);
                    fns{i} = LE{4};
                    
                elseif strcmp(LS{i},CLL{5}) %L3
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = eigval(:,:,:,3);
                    fns{i} = LE{5};
                    
                elseif strcmp(LS{i},CLL{6}) %RD
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = (eigval(:,:,:,2) + eigval(:,:,:,3))/2;
                    fns{i} = LE{6};
                    
                elseif strcmp(LS{i},CLL{7}) %RA
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    data_v{i} = EDTI_Library.E_DTI_RA(eigval);
                    fns{i} = LE{7};
                    
                elseif strcmp(LS{i},CLL{8}) %B0
                    
                    load(fin,'DWIB0')
                    
                    if ~exist('DWIB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = DWIB0;
                    fns{i} = LE{8};
                    
                elseif strcmp(LS{i},CLL{9}) %DT
                    
                    load(fin,'DT')
                    
                    if ~exist('DT','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DT = EDTI_Library.E_DTI_DT_mat2cell(DT);
                    data_v{i} = cat(4,DT{1},DT{2},DT{3},DT{4},DT{5},DT{6});
                    clear DT;
                    fns{i} = LE{9};
                    
                elseif strcmp(LS{i},CLL{10}) %CL
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(eigval);
                    data_v{i} = CL;
                    fns{i} = LE{10};
                    
                elseif strcmp(LS{i},CLL{11}) %CP
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(eigval);
                    data_v{i} = CP;
                    fns{i} = LE{11};
                    
                elseif strcmp(LS{i},CLL{12}) %CS
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(eigval);
                    data_v{i} = CS;
                    fns{i} = LE{12};
                    
                elseif strcmp(LS{i},CLL{13}) %FE
                    
                    load(fin,'FE')
                    
                    if ~exist('FE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = FE;
                    fns{i} = LE{13};
                    
                elseif strcmp(LS{i},CLL{14}) %SE
                    
                    load(fin,'SE')
                    
                    if ~exist('SE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = SE;
                    fns{i} = LE{14};
                    
                elseif strcmp(LS{i},CLL{15}) %TE
                    
                    load(fin,'FE','SE')
                    
                    if ~exist('FE','var') || ~exist('SE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    TE = FE;
                    
                    TE(:,:,:,1) = FE(:,:,:,2).*SE(:,:,:,3);
                    TE(:,:,:,1) = TE(:,:,:,1) - FE(:,:,:,3).*SE(:,:,:,2);
                    TE(:,:,:,2) = FE(:,:,:,3).*SE(:,:,:,1);
                    TE(:,:,:,2) = TE(:,:,:,2) - FE(:,:,:,1).*SE(:,:,:,3);
                    TE(:,:,:,3) = FE(:,:,:,1).*SE(:,:,:,2);
                    TE(:,:,:,3) = TE(:,:,:,3) - FE(:,:,:,2).*SE(:,:,:,1);
                    
                    data_v{i} = TE;
                    fns{i} = LE{15};
                    
                elseif strcmp(LS{i},CLL{16}) %abs(FE)
                    
                    load(fin,'FE')
                    
                    if ~exist('FE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    FE = FE(:,:,:,[2 1 3]);
                    data_v{i} = abs(FE);
                    fns{i} = LE{16};
                    
                elseif strcmp(LS{i},CLL{17}) %abs(SE)
                    
                    load(fin,'SE')
                    
                    if ~exist('SE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    SE = SE(:,:,:,[2 1 3]);
                    data_v{i} = abs(SE);
                    fns{i} = LE{17};
                    
                elseif strcmp(LS{i},CLL{18}) %abs(TE)
                    
                    load(fin,'FE','SE')
                    
                    if ~exist('FE','var') || ~exist('SE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    TE = FE;
                    
                    TE(:,:,:,1) = FE(:,:,:,2).*SE(:,:,:,3);
                    TE(:,:,:,1) = TE(:,:,:,1) - FE(:,:,:,3).*SE(:,:,:,2);
                    TE(:,:,:,2) = FE(:,:,:,3).*SE(:,:,:,1);
                    TE(:,:,:,2) = TE(:,:,:,2) - FE(:,:,:,1).*SE(:,:,:,3);
                    TE(:,:,:,3) = FE(:,:,:,1).*SE(:,:,:,2);
                    TE(:,:,:,3) = TE(:,:,:,3) - FE(:,:,:,2).*SE(:,:,:,1);
                    TE = TE(:,:,:,[2 1 3]);
                    data_v{i} = abs(TE);
                    fns{i} = LE{18};
                    
                elseif strcmp(LS{i},CLL{19}) %abs(FE)*FA
                    
                    load(fin,'FE','eigval')
                    
                    if ~exist('FE','var') || ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    
                    FA = EDTI_Library.FrAn_calc(eigval(:,:,:,1),eigval(:,:,:,2),eigval(:,:,:,3));
                    FE = FE(:,:,:,[2 1 3]);
                    data_v{i} = abs(FE).*(repmat(FA,[1 1 1 3]));
                    fns{i} = LE{19};
                    
                elseif strcmp(LS{i},CLL{20}) %abs(FE)*CL
                    
                    load(fin,'FE','eigval')
                    
                    if ~exist('FE','var') || ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    
                    [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(eigval);
                    FE = FE(:,:,:,[2 1 3]);
                    data_v{i} = abs(FE).*(repmat(CL,[1 1 1 3]));
                    fns{i} = LE{20};
                    
                elseif strcmp(LS{i},CLL{21}) %abs(TE)*CP
                    
                    load(fin,'FE','SE','eigval')
                    
                    if ~exist('FE','var') || ~exist('eigval','var') || ~exist('SE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    
                    [CL, CP, CS] = EDTI_Library.E_DTI_Westin_measures(eigval);
                    
                    TE = FE;
                    
                    TE(:,:,:,1) = FE(:,:,:,2).*SE(:,:,:,3);
                    TE(:,:,:,1) = TE(:,:,:,1) - FE(:,:,:,3).*SE(:,:,:,2);
                    TE(:,:,:,2) = FE(:,:,:,3).*SE(:,:,:,1);
                    TE(:,:,:,2) = TE(:,:,:,2) - FE(:,:,:,1).*SE(:,:,:,3);
                    TE(:,:,:,3) = FE(:,:,:,1).*SE(:,:,:,2);
                    TE(:,:,:,3) = TE(:,:,:,3) - FE(:,:,:,2).*SE(:,:,:,1);
                    TE = TE(:,:,:,[2 1 3]);
                    data_v{i} = abs(TE).*(repmat(CP,[1 1 1 3]));
                    fns{i} = LE{21};
                    
                elseif strcmp(LS{i},CLL{22}) %mean_res_DWI
                    
                    load(fin,'DT','DWI','DWIB0','b','NrB0','MDims')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DT','var') || ~exist('b','var') || ~exist('MDims','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    clear DWI DT DWIB0;
                    
                    [m_res_DWI, m_res_B0, M_res_DWI, M_res_B0] = EDTI_Library.E_DTI_Get_residual_info(fin);
                    
                    data_v{i} = m_res_DWI;
                    fns{i} = LE{22};
                    
                elseif strcmp(LS{i},CLL{23}) %max_res_DWI
                    
                    load(fin,'DT','DWI','DWIB0','b','NrB0','MDims')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DT','var') || ~exist('b','var') || ~exist('MDims','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    clear DWI DT DWIB0;
                    
                    [m_res_DWI, m_res_B0, M_res_DWI, M_res_B0] = EDTI_Library.E_DTI_Get_residual_info(fin);
                    
                    data_v{i} = M_res_DWI;
                    fns{i} = LE{23};
                    
                elseif strcmp(LS{i},CLL{24}) %mean_res_B0
                    
                    load(fin,'DT','DWI','DWIB0','b','NrB0','MDims')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DT','var') || ~exist('b','var') || ~exist('MDims','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    clear DWI DT DWIB0;
                    
                    [m_res_DWI, m_res_B0, M_res_DWI, M_res_B0] = EDTI_Library.E_DTI_Get_residual_info(fin);
                    
                    data_v{i} = m_res_B0;
                    fns{i} = LE{24};
                    
                elseif strcmp(LS{i},CLL{25}) %max_res_B0
                    
                    load(fin,'DT','DWI','DWIB0','b','NrB0','MDims')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DT','var') || ~exist('b','var') || ~exist('MDims','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    clear DWI DT DWIB0;
                    
                    [m_res_DWI, m_res_B0, M_res_DWI, M_res_B0] = EDTI_Library.E_DTI_Get_residual_info(fin);
                    
                    data_v{i} = M_res_B0;
                    fns{i} = LE{25};
                    
                elseif strcmp(LS{i},CLL{26}) %std_DWI
                    
                    load(fin,'DWI','NrB0')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_std_DWI(DWI, NrB0);
                    clear DWI;
                    fns{i} = LE{26};
                    
                elseif strcmp(LS{i},CLL{27}) %std_B0
                    
                    load(fin,'DWI','NrB0')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_std_B0(DWI, NrB0);
                    clear DWI;
                    fns{i} = LE{27};
                    
                elseif strcmp(LS{i},CLL{28}) %skew_L
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    mask = ~isnan(eigval(:,:,:,1));
                    data_v{i} = EDTI_Library.E_DTI_skewness_eigval(eigval,mask);
                    fns{i} = LE{28};
                    
                elseif strcmp(LS{i},CLL{29}) %PNPIV
                    
                    load(fin,'DWI','NrB0','eigval')
                    mask_PIS = ~isnan(eigval(:,:,:,1));
                    if ~exist('DWI','var') || ~exist('NrB0','var') || ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    [M_, M1_, M2_] = EDTI_Library.E_DTI_phpv(DWI, NrB0);
                    clear DWI;
                    
                    eigval = abs(eigval);
                    eigval = sort(eigval,4,'descend');
                    
                    FA = FrAn_calc(eigval(:,:,:,1),eigval(:,:,:,2),eigval(:,:,:,3));
                    clear eigval;
                    
                    FA_ = FA;
                    FA_(M_>0)=0;
                    FA(M_>0)=1;
                    
                    FA(~mask_PIS)=nan;
                    FA_(~mask_PIS)=nan;
                    
                    data_v{i} = cat(4,FA,FA_,FA_);
                    fns{i} = LE{29};
                    
                elseif strcmp(LS{i},CLL{30}) %DWIs
                    
                    load(fin,'DWI')
                    
                    if ~exist('DWI','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    %         for pp=1:length(DWI)
                    %            DWI{pp} = DWI{pp}(:,:,36:38);
                    %         end
                    
                    
                    data_v{i} = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
                    clear DWI;
                    fns{i} = LE{30};
                    
                elseif strcmp(LS{i},CLL{31}) %Kappa
                    
                    load(fin,'FE')
                    
                    if ~exist('FE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    [Kappa, IVDC] = EDTI_Library.E_DTI_IVDC(FE,1);
                    data_v{i} = Kappa;
                    clear FE;
                    fns{i} = LE{31};
                    
                elseif strcmp(LS{i},CLL{32}) %IVDC
                    
                    load(fin,'FE')
                    
                    if ~exist('FE','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    [Kappa, IVDC] = EDTI_Library.E_DTI_IVDC(FE,1);
                    data_v{i} = IVDC;
                    clear FE;
                    fns{i} = LE{32};
                    
                elseif strcmp(LS{i},CLL{33}) %B-matrix
                    
                    load(fin,'b')
                    
                    if ~exist('b','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    data_v{i} = b;
                    clear b;
                    fns{i} = LE{33};
                    
                elseif strcmp(LS{i},CLL{34}) %gradient directions
                    
                    load(fin,'g')
                    
                    if ~exist('g','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    data_v{i} = g;
                    clear g;
                    fns{i} = LE{34};
                    
                elseif strcmp(LS{i},CLL{35}) %ADCs
                    
                    load(fin,'bval','DWI','NrB0','DWIB0')
                    
                    if ~exist('DWI','var') || ~exist('bval','var') || ~exist('NrB0','var') || ~exist('DWIB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI = EDTI_Library.E_DTI_DWI_2_ADC(DWI,bval,NrB0,DWIB0); % is not really 'DWI' (memory efficiency ;-)
                    
                    data_v{i} = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
                    
                    clear DWI;
                    fns{i} = LE{35};
                    
                elseif strcmp(LS{i},CLL{36}) %DT_ADCs
                    
                    load(fin,'bval','DT','g')
                    
                    if ~exist('DT','var') || ~exist('bval','var') || ~exist('g','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    ADC = EDTI_Library.E_DTI_DT_2_ADC(DT, g, bval);
                    ADC = EDTI_Library.E_DTI_DWI_cell2mat(ADC);
                    data_v{i} = ADC;
                    
                    clear ADC;
                    fns{i} = LE{36};
                    
                elseif strcmp(LS{i},CLL{37}) %DT_DWIs
                    
                    load(fin,'b','DT','NrB0','DWIB0')
                    
                    if ~exist('DT','var') || ~exist('b','var') || ~exist('NrB0','var') || ~exist('DWIB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI = EDTI_Library.E_DTI_DT_2_DWI_big(DT, b, DWIB0, NrB0);
                    clear DT DWIB0;
                    DWI = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
                    
                    data_v{i} = DWI;
                    
                    clear DWI;
                    fns{i} = LE{37};
                    
                elseif strcmp(LS{i},CLL{38}) %res_DWIs
                    
                    load(fin,'DT','NrB0','DWIB0','b','DWI')
                    
                    if ~exist('DT','var') || ~exist('b','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DWI','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI = EDTI_Library.E_DTI_residuals_DT(DT, b, DWIB0, NrB0, DWI); % is not really 'DWI' (memory efficiency ;-)
                    clear DT DWIB0;
                    DWI = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
                    
                    data_v{i} = DWI;
                    
                    clear DWI;
                    fns{i} = LE{38};
                    
                elseif strcmp(LS{i},CLL{39}) %res_ADCs
                    
                    load(fin,'DT','NrB0','DWIB0','bval','DWI','g')
                    
                    if ~exist('DT','var') || ~exist('bval','var') || ~exist('NrB0','var') || ~exist('DWIB0','var') || ~exist('DWI','var') || ~exist('g','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI = EDTI_Library.E_DTI_ADC_res_DT(DT, g, bval, DWI, NrB0, DWIB0); % is not really 'DWI' (memory efficiency ;-)
                    clear DT DWIB0;
                    DWI = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
                    
                    data_v{i} = DWI;
                    
                    clear DWI;
                    fns{i} = LE{39};
                    
                elseif strcmp(LS{i},CLL{40}) %mask
                    
                    load(fin,'FA')
                    
                    if ~exist('FA','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = single(~isnan(FA));
                    fns{i} = LE{40};
                    
                elseif strcmp(LS{i},CLL{41}) %Lattice Index
                    
                    load(fin,'DT')
                    
                    if ~exist('DT','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DT = EDTI_Library.E_DTI_DT_mat2cell(DT);
                    data_v{i} = EDTI_Library.E_DTI_LI(DT, 1);
                    fns{i} = LE{41};
                    clear DT;
                    
                elseif strcmp(LS{i},CLL{42}) %Mean DWIs
                    
                    load(fin,'DWI','NrB0')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI = EDTI_Library.E_DTI_mean_DWI(DWI, NrB0);
                    data_v{i} = DWI;
                    fns{i} = LE{42};
                    clear DWI;
                    
                elseif strcmp(LS{i},CLL{43}) %Mean B0s
                    
                    load(fin,'DWI','NrB0')
                    
                    if ~exist('DWI','var') || ~exist('NrB0','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    DWI(NrB0+1:end)=[];
                    DWI = EDTI_Library.E_DTI_mean_DWI(DWI, 0);
                    data_v{i} = DWI;
                    fns{i} = LE{43};
                    clear DWI;
                    
                elseif strcmp(LS{i},CLL{44}) %HARDI
                    
                    dummy = EDTI_Library.E_DTI_Get_HARDI_stuff(fin);
                    if isempty(dummy)
                        continue;
                    end
                    data_v{i} = dummy;
                    fns{i} = LE{44};
                    
                elseif strcmp(LS{i},CLL{45}) %KT
                    
                    load(fin,'KT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    KT = EDTI_Library.E_DTI_DWI_cell2mat(KT);
                    
                    data_v{i} = KT;
                    fns{i} = LE{45};
                    
                elseif strcmp(LS{i},CLL{46}) %MK
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    load(fin,'info')
                    if isfield(info,'MK')
                        data_v{i} = info.MK;
                    else
                        data_v{i} = EDTI_Library.E_DTI_Mean_Kurtosis_c(KT,DT);
                    end
                    fns{i} = LE{46};
                    
                elseif strcmp(LS{i},CLL{47}) %AK
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    load(fin,'info')
                    if isfield(info,'K_para')
                        data_v{i} = info.K_para;
                    else
                        data_v{i} = EDTI_Library.E_DTI_Axial_Kurtosis_c(KT,DT);
                    end
                    fns{i} = LE{47};
                    
                elseif strcmp(LS{i},CLL{48}) %RK
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    load(fin,'info')
                    if isfield(info,'K_perp')
                        data_v{i} = info.K_perp;
                    else
                        data_v{i} = EDTI_Library.E_DTI_Radial_Kurtosis_c(KT,DT);
                    end
                    fns{i} = LE{48};
                    
                elseif strcmp(LS{i},CLL{49}) %KA
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    load(fin,'info')
                    if isfield(info,'K_an')
                        data_v{i} = info.K_an;
                    else
                        data_v{i} = EDTI_Library.E_DTI_Kurtosis_Anisotropy_c(KT,DT);
                    end
                    fns{i} = LE{49};
                    
                elseif strcmp(LS{i},CLL{50}) %AM
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_Anisotropy_Mode(eigval);
                    fns{i} = LE{50};
                    
                elseif strcmp(LS{i},CLL{51}) %DTN
                    
                    load(fin,'eigval')
                    
                    if ~exist('eigval','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_Tensor_Norm(eigval);
                    fns{i} = LE{51};
                    
                elseif strcmp(LS{i},CLL{52}) %FA Skeleton
                    
                    load(fin,'FA','VDims')
                    
                    if ~exist('FA','var') || ~exist('VDims','var')
                        disp(['Format of ''' fin ''' not correct, skipping data!'])
                        continue;
                    end
                    
                    FA(isnan(FA))=0;
                    FA=FA/sqrt(3);
                    FA(FA>1)=1;
                    FA(FA<0)=0;
                    
                    data_v{i} = EDTI_Library.E_DTI_FA_Skeletonization(FA,VDims);
                    fns{i} = LE{52};
                    
                elseif strcmp(LS{i},CLL{53}) %AWF
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_AWF_DKI(KT,DT);
                    
                    
                    fns{i} = LE{53};
                    
                elseif strcmp(LS{i},CLL{54}) %TORT
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_EAS_DKI(KT,DT);
                    
                    fns{i} = LE{54};
                    
                elseif strcmp(LS{i},CLL{55}) %AxEAD
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_axial_EAS_DKI(KT,DT);
                    
                    fns{i} = LE{55};
                    
                elseif strcmp(LS{i},CLL{56}) %RadEAD
                    
                    load(fin,'KT','DT')
                    
                    if ~exist('KT','var')
                        disp(['Format of ''' fin ''' not correct for DKI stuff, skipping data!'])
                        continue;
                    end
                    
                    data_v{i} = EDTI_Library.E_DTI_radial_EAS_DKI(KT,DT);
                    
                    fns{i} = LE{56};
                    
                end
                
            end
            
            warning on all
            
            for j=1:length(LS)
                if ~isempty(data_v{j})
                    if strcmp(fns{j}(end-2:end),'nii')
                        data_v{j}(isnan(data_v{j}))=0;
                        EDTI_Library.E_DTI_write_nifti_file(data_v{j},VDims,[folout filesep NAME fns{j}])
                        data_v{j}=[];
                    else
                        EDTI_Library.E_DTI_Write_text(data_v{j},[folout filesep NAME fns{j}])
                    end
                end
            end
            
            disp('Done!')
            
        end
        
        % From ExploreDTI: Generates seed points for tractography in a given a mask
        function seedPoint = E_DTI_Get_Seeds_WBT(seed_mask, SR, VDims, params)
            
            S = size(seed_mask);
            x = single(VDims(1):VDims(1):VDims(1)*S(1));
            lx = length(x);
            y = single(VDims(2):VDims(2):VDims(2)*S(2));
            ly = length(y);
            z = single(VDims(3):VDims(3):VDims(3)*S(3));
            lz = length(z);
            
            X = repmat(x',[1 ly lz]);
            Y = repmat(y,[lx 1 lz]);
            p = zeros(1,1,lz);
            p(:,:,:) = z;
            Z = repmat(p,[lx ly 1]);
            clear x y z p lx ly lz;
            
            x_min = min(X(seed_mask));
            x_max = max(X(seed_mask));
            y_min = min(Y(seed_mask));
            y_max = max(Y(seed_mask));
            z_min = min(Z(seed_mask));
            z_max = max(Z(seed_mask));
            clear X Y Z
            
            if SR(1)>x_max || SR(2)>y_max || SR(3)>z_max
                seedPoint = [];
                return;
            end
            
            
            xrange = x_min:SR(1):x_max;
            yrange = y_min:SR(2):y_max;
            zrange = z_min:SR(3):z_max;
            
            seedPoint = zeros(3,size(xrange,2)*size(yrange,2)*size(zrange,2), 'single');
            i = 0;
            for x = xrange
                for y = yrange
                    for z = zrange
                        
                        if params.randp==1
                            xx = x + SR(1)*(rand-0.5);
                            yy = y + SR(2)*(rand-0.5);
                            zz = z + SR(3)*(rand-0.5);
                        else
                            xx = x;
                            yy = y;
                            zz = z;
                        end
                        
                        v = floor([xx yy zz]./VDims);
                        if all([v>0 v<=size(seed_mask)])
                            if(seed_mask(v(1),v(2),v(3)))
                                i = i + 1;
                                seedPoint(:,i) = [xx yy zz]';
                            end
                        end
                    end
                end
            end
            seedPoint = seedPoint(:,1:i);
        end
        
        % This function is used to track multiple FODs (mFOD). Written by A. De
        % Luca
        function WholeBrainTracking_mDRL_fast_exe(f_in, suffix, f_out, p, use_linear_or_majority)
            
            if(strcmp(use_linear_or_majority,'linear'))
                use_linear_or_majority = 1;
            elseif(strcmp(use_linear_or_majority,'majority'))
                use_linear_or_majority = 2;
            else
                error('Unexpected interpolation.');
            end
            
            disp('Performing mFOD tracking (FOD interpolation) for file:')
            disp(['' f_in ''])
            
            disp('Calculating FOD...')
            
            [fractions,VD] = EDTI_Library.E_DTI_read_nifti_file([suffix '_fractions.nii']);
            [~,IX] = max(fractions(:,:,:,:),[],4);
            
            CSD_FOD = [];
            to_eliminate = zeros(size(fractions,4),1);
            
            for ij=1:size(fractions,4)
                try
                    VOL = EDTI_Library.E_DTI_read_nifti_file([suffix '_CSD_FOD_' num2str(ij) '.nii']);
                    if(isempty(CSD_FOD))
                        CSD_FOD = zeros([size(VOL) size(fractions,4)]);
                    end
                    if(isempty(VOL))
                        disp(['Found ' num2str(ij-1) ' FODs']);
                        CSD_FOD = CSD_FOD(:,:,:,:,1:ij-1);
                        break;
                    end
                    CSD_FOD(:,:,:,:,ij) = VOL;
                catch
                    CSD_FOD(:,:,:,:,ij) = 0;
                    to_eliminate(ij) = 1;
                end
            end
            
            % if(to_eliminate == size(fractions,4))
            %     % there is an isotropic compartment
            %     BadVoxels = fractions(:,:,:,end) > 0.9;%
            %     %     BadVoxels = IX == to_eliminate;
            %     for ij=1:size(CSD_FOD,5)
            %         for ih=1:size(CSD_FOD,4)
            %             V = CSD_FOD(:,:,:,ih,ij);
            %             V(BadVoxels) = 0;
            %             CSD_FOD(:,:,:,ih,ij) = V;
            %         end
            %     end
            % end
            
            fractions(:,:,:,to_eliminate == 1) = [];
            CSD_FOD(:,:,:,:,to_eliminate == 1) = [];
            
            [sx,sy,sz,sh,sf] = size(CSD_FOD);
            CSD_FOD = reshape(CSD_FOD,sx*sy*sz,sh,sf);
            [~,IX] = max(fractions(:,:,:,1:sf),[],4);
            
            if(use_linear_or_majority == 1)
                
                % Save a linterp FOD
                CSD_FOD_tmp = zeros([sx*sy*sz sh]);
                fr_col = reshape(fractions(:,:,:,1:sf),sx*sy*sz,sf);
                for ij=1:sf
                    the_norm = 1;%prctile(CSD_FOD(fr_col(:,ij)>0.7,1,1),95)/prctile(CSD_FOD(fr_col(:,1)>0.7,1,ij),95);
%                     the_norm = max(nanmean(CSD_FOD(IX==1,:,1)))/max(nanmean(CSD_FOD(IX==ij,:,ij)));
                    for ik=1:size(fr_col,1)
                        if(fr_col(ik,ij) == 0 || CSD_FOD(ik,1,ij) == 0 || ~isfinite(CSD_FOD(ik,1,ij)))
                            continue
                        end
                        CSD_FOD_tmp(ik,:) = CSD_FOD_tmp(ik,:) + fr_col(ik,ij)*CSD_FOD(ik,:,ij)*the_norm;
                    end
                end
                CSD_FOD_tmp = reshape(CSD_FOD_tmp,sx,sy,sz,sh);
                
                FOD_1 = CSD_FOD_tmp(:,:,:,1);
                FA_k = load(f_in,'FA');
                ratio = 1;%0.1/mean(FOD_1(FA_k.FA(:)>0.6));
                %     disp(['Ratio is: ' num2str(ratio)]);
                CSD_FOD_tmp = single(CSD_FOD_tmp*ratio);
                
                EDTI_Library.E_DTI_write_nifti_file(CSD_FOD_tmp,VD,[suffix '_CSD_FOD_linterp.nii']);
                
            elseif(use_linear_or_majority == 2)
                
                % Save a majority FOD
                CSD_FOD_tmp = zeros([sx*sy*sz sh]);
                fr_col = reshape(fractions(:,:,:,1:sf),sx*sy*sz,sf);
                the_norm = ones(sf,1);
                
                for ij=1:sf
                    the_norm(ij) = 1;%prctile(CSD_FOD(fr_col(:,ij)>0.7,1,1),95)/prctile(CSD_FOD(fr_col(:,1)>0.7,1,ij),95);
                end
                
                for ij=1:sf
                    SEL = IX == ij;
                    CSD_FOD_tmp(SEL,:) = CSD_FOD(SEL,:,ij)*the_norm(ij);
                end
                CSD_FOD_tmp = reshape(CSD_FOD_tmp,sx,sy,sz,sh);
                
                FOD_1 = CSD_FOD_tmp(:,:,:,1);
                FA_k = load(f_in,'FA');
                ratio = 1;%0.1/mean(FOD_1(FA_k.FA(:)>0.6));
                %     disp(['Ratio is: ' num2str(ratio)]);
                CSD_FOD_tmp = CSD_FOD_tmp*ratio;
                
                EDTI_Library.E_DTI_write_nifti_file(CSD_FOD_tmp,VD,[suffix '_CSD_FOD_majority.nii']);
                
            else
                disp('Unsupported');
                return;
            end
            
            EDTI_Library.WholeBrainFODTractography(f_in,CSD_FOD_tmp,p,f_out)
            
        end
        
        % Wraooer around the Gibbs ringing correction
        function GibbsRingingCorrection(fnp_in,fnp_out)
            
            try
                txt_file = load([fnp_in(1:end-4) '.txt']);
            catch
                error(['Cannot find a .txt gradient file associated to ' fnp_in]);
            end
            
            % bval = sum(txt_file(:,[1 4 6]),2);
            [~,~,IX] = EDTI.GetNumOfB0sDWIs('nii_file',fnp_in);
            % IX = find(bval < 1);
            IXd = diff(IX);
            if(~isempty(IXd) && any(IXd) ~= 1)
                error('The .nii file should be sorted per ascending b-value');
            end
            
            nrb0 = EDTI.GetNumOfB0sDWIs('nii_file',fnp_in);
            
            parameters.NrB0 = nrb0;
            parameters.lambda = 100;
            parameters.iter = 100;
            parameters.ss = 0.01;
            parameters.ip = 3;
            
            parameters.ext = '';
            
            
            EDTI_Library.E_DTI_Gibbs_Ringing_removal_with_TV_exe(fnp_in,fnp_out,parameters);
            
        end
        
        % From ExploreDTI: perform Gibbs ringing correction
        function suc  = E_DTI_Gibbs_Ringing_removal_with_TV_exe(f_in,f_out,parameters)
            
            suc = 1;
            
            try
                [DWI, VDims] = EDTI_Library.E_DTI_read_nifti_file(f_in);
            catch
                suc = 0;
                disp('Could not load file:')
                disp(f_in)
                return;
            end
            
            if ndims(DWI)~=4 && ndims(DWI)~=3
                suc = 0;
                disp('Error for file:')
                disp(f_in)
                disp('Number of dimensions should be 3 or 4!')
                return;
            end
            
            DWI = EDTI_Library.E_DTI_DWI_mat2cell(DWI);
            
            DWI(1:parameters.NrB0) = EDTI_Library.E_DTI_Gibbs_Ringing_removal_TV(DWI(1:parameters.NrB0),...
                parameters.lambda,parameters.iter,parameters.ss,parameters.ip);
            
            DWI = EDTI_Library.E_DTI_DWI_cell2mat(DWI);
            
            try
                EDTI_Library.E_DTI_write_nifti_file(DWI,VDims,f_out);
            catch
                suc = 0;
                disp('Error for file:')
                disp(f_in)
                disp('Probably memory issues...')
                return;
            end
            
        end
        
        % From ExploreDTI: The actual gibbs ringing correction
        function B0 = E_DTI_Gibbs_Ringing_removal_TV(B0,lambda,timeSteps,stepsize,recon)
            
            cla = class(B0{1});
            
            for i=1:length(B0)
                B0{i}=double(B0{i});
            end
            
            if recon==1
                for i=1:length(B0);
                    B0{i} = permute(B0{i},[3 2 1]);
                end
            elseif recon==2
                for i=1:length(B0);
                    B0{i} = permute(B0{i},[1 3 2]);
                end
            end
            
            parfor i=1:length(B0)
                for j=1:size(B0{i},3)
                    dummy = double(B0{i}(:,:,j));
                    M = max(dummy(:));
                    dummy = dummy/M;
                    dummy = EDTI_Library.digitalTotalVariationFilter_2d(dummy,lambda,timeSteps,stepsize);
                    B0{i}(:,:,j) = dummy*M;
                end
            end
            
            if recon==1
                for i=1:length(B0);
                    B0{i} = permute(B0{i},[3 2 1]);
                end
            elseif recon==2
                for i=1:length(B0);
                    B0{i} = permute(B0{i},[1 3 2]);
                end
            end
            
            for i=1:length(B0)
                if strcmp(cla,'int16')
                    B0{i} = int16(round(B0{i}));
                elseif strcmp(cla,'uint16')
                    B0{i} = uint16(round(B0{i}));
                elseif strcmp(cla,'single')
                    B0{i} = single(B0{i});
                end
            end
            
        end
        
        % From ExploreDTI: Total variation filter
        function v = digitalTotalVariationFilter_2d(u0,lambda,number_iter,stepsize)
            
            % DIGITALTOTALVARIATIONFILTER_2D
            %
            % References: (1) Digitized {PDE} Method for Data Restoration. Ch. 16, p. 751 to 771,
            %                 in Analytic-Computational Methods in Applied Mathematics (2000), G. Anastassiou editor
            %             (2) Digital Total Variation Filtering as Postprocessing for Pseudospectral Methods for Conservation Laws,
            %                 Numerical Algorithms, vol. 41, p. 17-33, 2006
            % Inputs
            %       u0[][]     physical space function values (not spectral coefficients)
            %      lambda      fitting parameter
            %   number_iter      number of time marching steps
            % Output
            %        v[][]     the postprocessed function values
            % Called by:
            %   1) postProcessDriver2d.m
            % Notes:
            %   uses time-marching (Euler's method) to advance the nonlinear restoration to a steady state
            %   uses a 4 point neighborhood
            % Last modified: October 17, 2007
            
            
            a = 1e-4;
            a = a^2;
            N = length(u0(:,1));
            M = length(u0(1,:));
            v = u0; U = u0;
            
            dt = stepsize;
            
            k = 1;
            s = zeros(N,M);
            
            while k <= number_iter
                
                i = 2:N-1;
                j = 2:M-1;
                
                s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
                i=1;
                s(1,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 +          0           + (U(i+1,j)-U(i,j)).^2 + a );
                i = 2:N-1; j=M;
                s(i,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
                i=N; j = 2:M-1;
                s(N,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 +          0           + a );
                i = 2:N-1; j=1;
                s(i,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
                i=1; j=1;
                s(1,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 +          0           + (U(i+1,j)-U(i,j)).^2 + a );
                i=N; j=M;
                s(N,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           + (U(i-1,j)-U(i,j)).^2 +          0           + a );
                i=N; j=1;
                s(N,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 +          0           + a );
                i=1; j=M;
                s(1,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           +          0           + (U(i+1,j)-U(i,j)).^2 + a );
                
                %        v = euler(U,k*dt,dt,@F);
                
                
                i=2:N-1;
                j=2:M-1;
                
                v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                    ( U(i-1,j) - U(i,j) ).*( 1 + s(i,j)./s(i-1,j) ) + ...
                    ( U(i,j+1) - U(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                    ( U(i,j-1) - U(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                    lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   );
                %          end
                %        end
                
                v(1:N,1) = u0(1:N,1);
                v(1:N,M) = u0(1:N,M);
                v(1,1:M) = u0(1,1:M);
                v(N,1:M) = u0(N,1:M);
                
                k = k+1;  U = v;
                
            end  % while
            
        end
        
        % From ExploreDTI: Helper function for the automatic determination of
        % gradient signs / flips
        function f = PermutationOptimFun(perm,sign,p,logdt,mask,points)
            
            %
            R = eye(3,3); R = R(perm,:); R = bsxfun(@times,R,sign');
            R = EDTI_Library.DiffusionTensorRotationMatrix(R);
            logdt = R*logdt;
            logdt = mex_unvec(logdt,mask);
            tracts = EDTI_Library.track(p,logdt,points);
            f = -mean(tracts);
        end
        
        % From ExploreDTI: Helper function for the automatic determination of
        % gradient signs / flips
        function fval = SingleTrial(famask,p,n,perm_list,sign_list,logdt,mask)
            fval = zeros(size(perm_list,1),size(sign_list,1));
            points = EDTI_Library.SampleSeedPoints(famask, p.vdims, n);
            for i = 1:size(perm_list,1)
                perm = perm_list(i,:);
                for j = 1:size(sign_list,1)
                    sign = sign_list(j,:);
                    fval(i,j) = EDTI_Library.PermutationOptimFun(perm,sign,p,logdt,mask,points);
                end
            end
        end
        
        % From ExploreDTI: Helper function for the automatic determination of
        % gradient signs / flips
        function [grad] =  E_DTI_grad(bmat)
            bmat = double(bmat);
            IND = bmat(:,[1 4 6])<0;
            bmat([IND(:,1) false([size(IND,1) 2]) IND(:,2) false([size(IND,1) 1]) IND(:,2)])=0;
            bvals = sum(bmat(:,[1 4 6]),2);
            BSign_x   = sign(sign(bmat(:,1:3)) + 0.0001);
            BSign_y   = sign(sign(bmat(:,[2 4 5])) + 0.0001);
            Bo = bmat(:,1)==0;
            BSign = BSign_x;
            BSign([Bo Bo Bo]) = BSign_y([Bo Bo Bo]);
            grad = BSign .* sqrt(bmat(:,[1 4 6]));
            grad_n = sqrt(sum(grad.*grad,2));
            grad = grad./[grad_n grad_n grad_n];
            grad(isnan(grad))=0;
            
            if ~isreal(grad)
                grad
                disp('Warning: incorrect values of the B-matrix were encountered!')
                grad = real(grad);
                grad_n = sqrt(sum(grad.*grad,2));
                grad = grad./[grad_n grad_n grad_n];
                grad(isnan(grad))=0;
            end
            
            grad = cat(2,grad,bvals);
        end
        
        % From ExploreDTI: Helper function for the automatic determination of
        % gradient signs / flips
        function [flips, perms, correct, consistent] = E_DTI_Check_for_flip_perm_grads(mat_file_name)
            
            try
                
                correct = 1;
                flips = [];
                perms = [];
                
                load(mat_file_name,'VDims','b','DT')
                
                v2w = diag(VDims);
                v2w(4,4) = 1;
                grad = EDTI_Library.E_DTI_grad(b);
                
                for dwin = 1:length(DT)
                    DT{dwin} = smooth3(DT{dwin},'gaussian');
                end
                
                mask = ~isnan(DT{1});
                dt = EDTI_Library.vec_cell(DT,mask);
                
                grad(grad(:,4)>0,1:3) = bsxfun(@rdivide,grad(grad(:,4)>0,1:3),sqrt(sum(grad(grad(:,4)>0,1:3).^2,2)));
                curv = double(2);
                p.v2w = v2w;
                p.vdims = abs(p.v2w*[1 1 1 0]'-p.v2w*[2 2 2 0]'); p.vdims = p.vdims(1:3)';
                p.stepSize = 0.5*min(p.vdims);
                p.maxAngle = (180/pi)*2*asin(p.stepSize/(2*curv));
                p.lengthRange = [p.stepSize*3 p.stepSize*100];
                p.v2w = diag([p.vdims 1]);
                
                [eigval, eigvec] = mex_dteig(double(dt));
                
                logeigval = log(eigval);
                logdt = zeros(6,size(eigvec,2),class(eigval));
                logdt(1,:) = sum(eigvec([1 4 7],:).*logeigval.*eigvec([1 4 7],:));
                logdt(2,:) = sum(eigvec([1 4 7],:).*logeigval.*eigvec([2 5 8],:));
                logdt(3,:) = sum(eigvec([1 4 7],:).*logeigval.*eigvec([3 6 9],:));
                logdt(4,:) = sum(eigvec([2 5 8],:).*logeigval.*eigvec([2 5 8],:));
                logdt(5,:) = sum(eigvec([2 5 8],:).*logeigval.*eigvec([3 6 9],:));
                logdt(6,:) = sum(eigvec([3 6 9],:).*logeigval.*eigvec([3 6 9],:));
                fa = DTI.fa(eigval);
                meanfa = nanmean(fa);
                p.threshold = meanfa;
                famask = unvec(fa,mask) > meanfa;
                p.mask = famask;
                
                n = 50;
                
                perm_list = VChooseKO(1:3,3);
                sign_list = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1];
                
                
                trials = 10;
                for l = 1:8
                    fval = nan(size(perm_list,1),size(sign_list,1),trials);
                    prev_idx = nan;
                    consistent = true;
                    for k = 1:trials
                        fval(:,:,k) = EDTI_Library.SingleTrial(famask,p,n,perm_list,sign_list,logdt,mask);
                        fval_ = fval(:,:,k); fval_ = fval_(:);
                        [~, idx] = min(fval_);
                        if isnan(prev_idx)
                            prev_idx = idx;
                        end
                        if idx ~= prev_idx
                            consistent = false;
                            break;
                        end
                    end
                    if consistent
                        break;
                    else
                        n = n * 2;
                    end
                end
                fval_ = nanmean(fval,3);
                [~,i] = min(fval_(:)); [a,b]=ind2sub(size(fval_),i);
                perms = perm_list(a,:);
                flips = sign_list(b,:);
                
                if ~all([perms flips]==[1 2 3 1 1 1])
                    correct=0;
                end
                
            catch me
                %     disp(me.message)
                flips = [];
                perms = [];
                correct = 1;
                consistent = [];
                return;
            end
        end
        
        % From ExploreDTI: Helper function for the resampling of 3D volumes
        function im = E_DTI_resample_nii_file(im, res, VDims)
            
            
            if ndims(im)==3
                
                try
                    im = EDTI_Library.E_DTI_resample_nii_file_3D(im, res, VDims);
                catch
                    im = [];
                    return;
                end
                
            elseif ndims(im)==4
                
                try
                    im1 = EDTI_Library.E_DTI_resample_nii_file_3D(im(:,:,:,1), res, VDims);
                catch
                    im = [];
                    return;
                end
                
                try
                    
                    im1 = repmat(im1,[1 1 1 size(im,4)]);
                    
                    for i=2:size(im,4);
                        
                        im1(:,:,:,i) = EDTI_Library.E_DTI_resample_nii_file_3D(im(:,:,:,i), res, VDims);
                        
                    end
                    
                    clear im;
                    im = im1;
                    clear im1;
                    
                catch
                    im = [];
                    return;
                end
            end
        end
        
        % From ExploreDTI: Helper function for the resampling of 3D volumes
        function im = E_DTI_resample_nii_file_3D(im, res, VDims)
            
            if isa(im,'single')
                fl = 1;
            elseif isa(im,'int16')
                fl = 2;
            elseif isa(im,'uint16')
                fl = 3;
            elseif isa(im,'double')
                fl = 4;
            elseif isa(im,'int8')
                fl = 5;
            elseif isa(im,'uint8')
                fl = 6;
            else
                fl = 7;
            end
            
            im = single(im);
            
            % res_method = 'linear';
            res_method = 'spline';
            
            recon_res = VDims;
            recon_mat = size(im);
            recon_dims = recon_mat.*recon_res;
            
            [xi,yi,zi] = ndgrid(recon_res(1):recon_res(1):recon_dims(1),...
                recon_res(2):recon_res(2):recon_dims(2),recon_res(3):recon_res(3):recon_dims(3));
            
            xi = single(xi);
            yi = single(yi);
            zi = single(zi);
            
            x_min = min(xi(:));
            x_max = max(xi(:));
            y_min = min(yi(:));
            y_max = max(yi(:));
            z_min = min(zi(:));
            z_max = max(zi(:));
            
            [xj,yj,zj] = ndgrid(x_min:res(1):x_max,...
                y_min:res(2):y_max,z_min:res(3):z_max);
            
            xj = single(xj);
            yj = single(yj);
            zj = single(zj);
            
            im = interpn(xi,yi,zi,im,xj,yj,zj,res_method);
            
            if fl==1
                im = single(im);
            elseif fl==2
                im = round(im);
                im(im>intmax('int16')) = intmax('int16');
                im = int16(im);
            elseif fl==3
                im = round(im);
                im(im>intmax('uint16')) = intmax('uint16');
                im = uint16(im);
            elseif fl==4
                im = double(im);
            elseif fl==5
                im = round(im);
                im(im>intmax('int8')) = intmax('int8');
                im = int8(im);
            elseif fl==6
                im = round(im);
                im(im>intmax('uint8')) = intmax('uint8');
                im = uint8(im);
            elseif fl==7
                im = single(im);
            end
        end
        
        function lmax = E_DTI_n2lmax(n)
            lmax = 2*(floor((sqrt(1+8*n)-3)/4));
        end
        
        function [y, mask] = vec_cell(x,mask)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            %
            if ~exist('mask','var')
                mask = ~isnan(x{1});
            end
            y = zeros([size(x,2) sum(mask(:))], class(x{1}));
            for k = 1:size(x,2);
                Dummy = x{k};
                y(k,:) = Dummy(mask(:));
            end
        end
        
        function y = unvec_cell(x,mask)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            %
            dims = [size(mask,1) size(mask,2) size(mask,3)];
            
            y = cell(1,size(x,1));
            
            for k = 1:size(x,1)
                if isfloat(x)
                    Dummy = NaN(dims, class(x));
                else
                    Dummy = zeros(dims, class(x));
                end
                Dummy(mask) = x(k,:);
                y{k} = Dummy;
            end
        end
        
        function seedPoints = SampleSeedPoints(mask, vdims, n)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            %
            seedPoints = zeros(3,n);
            i = 1;
            while i<=n
                point=(rand([3,1]).*(size(mask)-1)')+[1 1 1]';
                point_ = round(point);
                
                if mask(point_(1), point_(2), point_(3))==1
                    seedPoints(:,i)=point(:);
                    i=i+1;
                end
                
            end
            seedPoints = bsxfun(@times,seedPoints,vdims');
        end
        
        function M = DiffusionTensorRotationMatrix(R)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            M = [R(1,1)^2 2*R(1,1)*R(1,2) 2*R(1,1)*R(1,3) R(1,2)^2 2*R(1,2)*R(1,3) R(1,3)^2;
                R(1,1)*R(2,1) R(1,1)*R(2,2)+R(1,2)*R(2,1) R(1,1)*R(2,3)+R(1,3)*R(2,1) R(1,2)*R(2,2) R(1,2)*R(2,3)+R(1,3)*R(2,2) R(1,3)*R(2,3);
                R(1,1)*R(3,1) R(1,1)*R(3,2)+R(1,2)*R(3,1) R(1,1)*R(3,3)+R(1,3)*R(3,1) R(1,2)*R(3,2) R(1,2)*R(3,3)+R(1,3)*R(3,2) R(1,3)*R(3,3);
                R(2,1)^2 2*R(2,1)*R(2,2) 2*R(2,1)*R(2,3) R(2,2)^2 2*R(2,2)*R(2,3) R(2,3)^2;
                R(2,1)*R(3,1) R(2,1)*R(3,2)+R(2,2)*R(3,1) R(2,1)*R(3,3)+R(2,3)*R(3,1) R(2,2)*R(3,2) R(2,2)*R(3,3)+R(2,3)*R(3,2) R(2,3)*R(3,3);
                R(3,1)^2 2*R(3,1)*R(3,2) 2*R(3,1)*R(3,3) R(3,2)^2 2*R(3,2)*R(3,3) R(3,3)^2];
        end
        
        function tract = track(p,logdt,point)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            %
            point(4,:) = 1;
            voxel = p.v2w\point;
            voxel = voxel(1:3,:);
            point = point(1:3,:);
            logdt_m = mex_interp(logdt, voxel);
            
            mask = ~isnan(logdt_m(1,:));
            point = point(:,mask);
            logdt_m = logdt_m(:,mask);
            
            [logeigval, mydir] = mex_dteig(double(logdt_m));
            mydir = mydir(1:3,:);
            val = DTI.fa(exp(logeigval));
            
            mask = val > p.threshold;
            point = point(:,mask);
            mydir = mydir(:,mask);
            
            tract1 = EDTI_Library.trackOneDir(p,logdt,point, mydir);
            tract2 = EDTI_Library.trackOneDir(p,logdt,point,-mydir);
            
            tract = tract1+tract2;
            tract = tract.*p.stepSize;
        end
        
        function tract = trackOneDir(p,logdt,point,mydir)
            %
            % Copyright Ben Jeurissen (ben.jeurissen@uantwerpen.be)
            %
            %
            tract = zeros([1 size(point,2)]);
            flist = 1:size(point,2);
            
            for it = 1:(p.lengthRange(2)/p.stepSize)
                point = point + p.stepSize .* mydir;
                
                point(4,:) = 1;
                voxel = p.v2w\point;
                voxel = voxel(1:3,:);
                point = point(1:3,:);
                logdt_m = mex_interp(logdt, voxel);
                
                mask = ~isnan(logdt_m(1,:));
                point = point(:,mask);
                mydir = mydir(:,mask);
                logdt_m = logdt_m(:,mask);
                flist = flist(:,mask);
                
                [logeigval, newDir] = mex_dteig(double(logdt_m));
                newDir = newDir(1:3,:);
                eigval = exp(logeigval);
                val = DTI.fa(eigval);
                
                dot = sum(mydir.*newDir,1);
                myangle = (180/pi)*real(acos(abs(dot)));
                
                flipsign = sign(dot);
                mydir = flipsign([1 1 1],:).*newDir;
                
                mask = val > p.threshold | myangle < p.maxAngle;
                point = point(:,mask);
                mydir = mydir(:,mask);
                flist = flist(:,mask);
                val = val(mask);
                
                if isempty(point)
                    break
                end
                
                tract(flist) = tract(flist)+val;
            end
        end

    end
end
