%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function [yaw,roll,pitch,traslx,trasly,traslz,determinant] = DW_EddyCorrect_Win(dataname,correction_method,rotate_bvecs,reference_volume,show_plot,parameters_file)
if(nargin < 2)
    correction_method = 'eddy_correct';
end

if(exist('reference_volume','var') < 1)
    reference_volume = 1;
end

if(strcmp(correction_method,'eddy_correct'))
    system(['eddy_correct ' dataname ' ' dataname '_ec 0']);
elseif(strcmp(correction_method,'ants_affine'))
    
    if(exist([dataname '.nii.gz'],'file') > 0)
        data_ext = '.nii.gz';    
    else
        data_ext = '.nii';
    end
       
    dw_data = DW_LoadDataUnscale([dataname data_ext],[dataname '.bvec'],[dataname '.bval']);
    yaw = zeros(dw_data.hdr.dime.dim(5),1);
    pitch = zeros(dw_data.hdr.dime.dim(5),1);
    roll = zeros(dw_data.hdr.dime.dim(5),1);
    traslx = zeros(dw_data.hdr.dime.dim(5),1);
    trasly = zeros(dw_data.hdr.dime.dim(5),1);
    traslz = zeros(dw_data.hdr.dime.dim(5),1);
    determinant = zeros(dw_data.hdr.dime.dim(5),1);
    
    last_j = 0;
    
    ec_data.hdr = dw_data.hdr;
    ec_data.img = zeros(size(dw_data.img));
    ec_data.untouch = 1;
    ec_data.img(:,:,:,reference_volume) = dw_data.img(:,:,:,reference_volume);
    
    if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
        ec_data.bvecs = dw_data.bvecs;
    end
    
    ref_vol.hdr = dw_data.hdr;
    ref_vol.hdr.dime.dim(1) = 3;
    ref_vol.hdr.dime.dim(5) = 1;
    ref_vol.untouch = 1;
    if(exist('ANTS_TEMP','dir') < 1)
        mkdir('ANTS_TEMP');
        % SAVE REFERENCE VOLUME
        ref_vol.img = dw_data.img(:,:,:,reference_volume);
        save_untouch_nii(ref_vol,'ANTS_TEMP/ref_vol.nii.gz');
    else
        disp('ANTS_TEMP dir already exists, trying to recover last run');
        for j=1:length(dw_data.bvals)
            if(j == reference_volume)
                continue
            end
            try
                if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
                    % READ ROTATION MATRIX
                    f = fopen(sprintf('ANTS_TEMP/vol_%04dAffine.txt',j),'rt');
                    for waste=1:3
                        fgetl(f);
                    end
                    useful_line = fgetl(f);
                    fclose(f);
                    splitted_line = strsplit(useful_line,' ');
                    RMatrix = [str2num(splitted_line{2}) str2num(splitted_line{3}) str2num(splitted_line{4})
                        str2num(splitted_line{5}) str2num(splitted_line{6}) str2num(splitted_line{7})
                        str2num(splitted_line{8}) str2num(splitted_line{9}) str2num(splitted_line{10})];
                    ec_data.bvecs(j,:) = RMatrix*(dw_data.bvecs(j,:)');
                    ec_data.bvecs(j,:) = ec_data.bvecs(j,:)/norm(ec_data.bvecs(j,:),2);
                    yaw(j) = atan2d(RMatrix(2,1),RMatrix(1,1));
                    pitch(j) = atan2d(-RMatrix(3,1),sqrt(RMatrix(3,2).^2+RMatrix(3,3).^2));
                    roll(j) = atan2d(RMatrix(3,2),RMatrix(3,3));
                    traslx(j) = str2num(splitted_line{10});
                    trasly(j) = str2num(splitted_line{11});
                    traslz(j) = str2num(splitted_line{12});
                    determinant(j) = det(RMatrix);
                end
                temp_vol = load_untouch_nii(sprintf('ANTS_TEMP/warped_vol_%04d.nii.gz',j));
                ec_data.img(:,:,:,j) = temp_vol.img;
                last_j = j;
            catch
                break
            end
        end
    end
    
    
    
    if(exist('show_plot','var') > 0 && show_plot == 1)
        figure
    end
    % PROCESS ALL OTHER VOLUMES
    for j=last_j+1:length(dw_data.bvals)
        if(j == reference_volume)
            continue
        end
        ref_vol.img = dw_data.img(:,:,:,j);
        ref_vol.untouch = 1;
        save_untouch_nii(ref_vol,sprintf('ANTS_TEMP/vol_%04d.nii.gz',j));
        system(['ANTS 3 -m CC[ANTS_TEMP/ref_vol.nii.gz,' sprintf('ANTS_TEMP/vol_%04d.nii.gz',j) ',1,4] --number-of-affine-iterations 10000x1000x1000 -i 0 --affine-metric-type CC -o ' sprintf('ANTS_TEMP/vol_%04d',j)]);
        system(['WarpImageMultiTransform 3 ' sprintf('ANTS_TEMP/vol_%04d.nii.gz',j) ' ' sprintf('ANTS_TEMP/warped_vol_%04d.nii.gz',j) ' -R ANTS_TEMP/ref_vol.nii.gz ' sprintf('ANTS_TEMP/vol_%04dAffine.txt',j)]);
        temp_vol = load_untouch_nii(sprintf('ANTS_TEMP/warped_vol_%04d.nii.gz',j));
        ec_data.img(:,:,:,j) = temp_vol.img;
        
        if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
            % READ ROTATION MATRIX
            f = fopen(sprintf('ANTS_TEMP/vol_%04dAffine.txt',j),'rt');
            for waste=1:3
                fgetl(f);
            end
            useful_line = fgetl(f);
            fclose(f);
            splitted_line = strsplit(useful_line,' ');
            RMatrix = [str2num(splitted_line{2}) str2num(splitted_line{3}) str2num(splitted_line{4})
                str2num(splitted_line{5}) str2num(splitted_line{6}) str2num(splitted_line{7})
                str2num(splitted_line{8}) str2num(splitted_line{9}) str2num(splitted_line{10})];
            ec_data.bvecs(j,:) = RMatrix*(dw_data.bvecs(j,:)');
            ec_data.bvecs(j,:) = ec_data.bvecs(j,:)/norm(ec_data.bvecs(j,:),2);
            yaw(j) = atan2d(RMatrix(2,1),RMatrix(1,1));
            pitch(j) = atan2d(-RMatrix(3,1),sqrt(RMatrix(3,2).^2+RMatrix(3,3).^2));
            roll(j) = atan2d(RMatrix(3,2),RMatrix(3,3));
            traslx(j) = str2num(splitted_line{10});
            trasly(j) = str2num(splitted_line{11});
            traslz(j) = str2num(splitted_line{12});
            determinant(j) = det(RMatrix);
        end
        
        if(exist('show_plot','var') > 0 && show_plot == 1)
            subplot(211)
            hold off
            plot(yaw,'r');
            hold on
            plot(pitch,'g');
            plot(roll,'b');
            title('Yaw pitch roll');
            ylabel('deg');
            legend({'RotX','RotY','RotZ'});
            subplot(212)
            hold off
            plot(traslx,'r');
            hold on
            plot(trasly,'g');
            plot(traslz,'b');
            ylabel('mm');
            title('Traslations X-Y-Z');
            legend({'X','Y','Z'});
            drawnow
        end
    end
    save_untouch_nii(ec_data,[dataname '_ants_ec.nii.gz']);
    if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
        bvecs = ec_data.bvecs';
        save([dataname '_ants_ec.bvec'],'bvecs', '-ascii');
    end
    
    if(exist('show_plot','var') > 0 && show_plot == 1)
        savefig([dataname '_ants.fig']);
    end
    save([dataname '_ants'],'yaw','roll','pitch','determinant','traslx','trasly','traslz');
    rmdir('ANTS_TEMP','s');
elseif(strcmp(correction_method,'elastix_affine'))
%     elastix_path = '/data/home/adeluca/elastix_build/';
    [~,ecmd] = system('which elastix');
    ecmd = ecmd(1:end-1);
    ecmd_parts = strsplit(ecmd,'/');
    ecmd_path = '/';
    for ijk=1:length(ecmd_parts)-1
        ecmd_path = fullfile(ecmd_path,ecmd_parts{ijk});
    end
    elastix_path=ecmd_path;%'/home/adeluca/Tools/elastix/bin/';
%     elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' /data/home/adeluca/elastix_build/elastix'];
    elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' ecmd];% ' /home/adeluca/Tools/elastix/bin/elastix'];
    
    elastix_cmd = 'elastix_64';
    
    if(exist('parameters_file','var') > 0 && ~isempty(parameters_file))
        transform_file = which(parameters_file);
    else
        transform_file = which('elastix_affine_dti_linearinterp.txt');
    end
    if(isempty(transform_file))
        disp('Cannot find transformation file');
        return
    end
    if(exist([dataname '.nii.gz'],'file') > 0)
           data_ext = '.nii.gz';    
       else
           data_ext = '.nii';
    end
    dw_data = DW_LoadDataUnscale([dataname data_ext],[dataname '.bvec'],[dataname '.bval']);
    yaw = zeros(dw_data.hdr.dime.dim(5),1);
    pitch = zeros(dw_data.hdr.dime.dim(5),1);
    roll = zeros(dw_data.hdr.dime.dim(5),1);
    traslx = zeros(dw_data.hdr.dime.dim(5),1);
    trasly = zeros(dw_data.hdr.dime.dim(5),1);
    traslz = zeros(dw_data.hdr.dime.dim(5),1);
    determinant = zeros(dw_data.hdr.dime.dim(5),1);
    
    last_j = 0;
    
    ec_data.hdr = dw_data.hdr;
    ec_data.img = zeros(size(dw_data.img));
    ec_data.untouch = 1;
    ec_data.img(:,:,:,reference_volume) = dw_data.img(:,:,:,reference_volume);
    
    if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
        ec_data.bvecs = dw_data.bvecs;
    end
    
    ref_vol.hdr = dw_data.hdr;
    ref_vol.hdr.dime.dim(1) = 3;
    ref_vol.hdr.dime.dim(5) = 1;
    ref_vol.untouch = 1;
    if(exist('ELASTIX_TEMP','dir') < 1)
        mkdir('ELASTIX_TEMP');
        % SAVE REFERENCE VOLUME
        ref_vol.img = dw_data.img(:,:,:,reference_volume);
        save_untouch_nii(ref_vol,'ELASTIX_TEMP/ref_vol.nii.gz');
    else
        disp('ELASTIX_TEMP dir already exists, trying to recover last run');
        for j=1:length(dw_data.bvals)
            if(j == reference_volume)
                continue
            end
            try
                if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
                    % READ ROTATION MATRIX
                    f = fopen(sprintf('ELASTIX_TEMP/vol_%04d/TransformParameters.0.txt',j),'rt');
                    while(isempty(strfind(fgetl(f),'CenterOfRotationPoint')))
                        if(feof(f))
                            disp('Error reading transformation file');
                            fclose(f);
                            return
                        end
                    end
                    useful_line = fgetl(f); % transform matrix
                    useful_line = strrep(useful_line,'(','');
                    useful_line = strrep(useful_line,')','');
                    fclose(f);
                    splitted_line = strsplit(useful_line,' ');
                    RMatrix = [str2num(splitted_line{2}) str2num(splitted_line{3}) str2num(splitted_line{4})
                        str2num(splitted_line{5}) str2num(splitted_line{6}) str2num(splitted_line{7})
                        str2num(splitted_line{8}) str2num(splitted_line{9}) str2num(splitted_line{10})];
                    ec_data.bvecs(j,:) = RMatrix*(dw_data.bvecs(j,:)');
                    ec_data.bvecs(j,:) = ec_data.bvecs(j,:)/norm(ec_data.bvecs(j,:),2);
                    yaw(j) = atan2d(RMatrix(2,1),RMatrix(1,1));
                    pitch(j) = atan2d(-RMatrix(3,1),sqrt(RMatrix(3,2).^2+RMatrix(3,3).^2));
                    roll(j) = atan2d(RMatrix(3,2),RMatrix(3,3));
                    traslx(j) = str2num(splitted_line{11});
                    trasly(j) = str2num(splitted_line{12});
                    traslz(j) = str2num(splitted_line{13});
                    determinant(j) = det(RMatrix);
                end
                temp_vol = load_untouch_nii(sprintf('ELASTIX_TEMP/vol_%04d/result.0.nii',j));
                ec_data.img(:,:,:,j) = temp_vol.img;
                last_j = j;
            catch
                break
            end
        end
    end
    
    if(exist('show_plot','var') > 0 && show_plot == 1)
        figure
    end
    % PROCESS ALL OTHER VOLUMES
    for j=last_j+1:length(dw_data.bvals)
        if(j == reference_volume)
            temp_vol = load_untouch_nii('ELASTIX_TEMP/ref_vol.nii.gz');
            if(ec_data.hdr.dime.scl_slope ~= 0)
                temp_vol.img = temp_vol.img*ec_data.hdr.dime.scl_slope+ec_data.hdr.dime.scl_inter;
            end
            ec_data.hdr.dime.scl_slope = 1;
            ec_data.hdr.dime.scl_inter = 0;
            ec_data.img(:,:,:,j) = temp_vol.img;
            continue
        else
            ref_vol.img = dw_data.img(:,:,:,j);
            ref_vol.untouch = 1;
            save_untouch_nii(ref_vol,sprintf('ELASTIX_TEMP/vol_%04d.nii.gz',j));
            mkdir(sprintf('ELASTIX_TEMP/vol_%04d',j));
            system([elastix_cmd ' -f ELASTIX_TEMP/ref_vol.nii.gz -m ' sprintf('ELASTIX_TEMP/vol_%04d.nii.gz',j) ' -out ' sprintf('ELASTIX_TEMP/vol_%04d',j) ' -p ' transform_file]);
            %         system(['transformix -in ' sprintf('ELASTIX_TEMP/vol_%04d.nii.gz',j) ' -out ' sprintf('ELASTIX_TEMP/vol_%04d',j) ' -tp ' sprintf('ELASTIX_TEMP/vol_%04d/TransformParameters.0.txt',j)]);
            %         system(['ANTS 3 -m CC[TEMP/ref_vol.nii.gz,' sprintf('TEMP/vol_%04d.nii.gz',j) ',1,4] --number-of-affine-iterations 10000x1000x1000 -i 0 --affine-metric-type CC -o ' sprintf('TEMP/vol_%04d',j)]);
            %         system(['WarpImageMultiTransform 3 ' sprintf('TEMP/vol_%04d.nii.gz',j) ' ' sprintf('TEMP/warped_vol_%04d.nii.gz',j) ' -R TEMP/ref_vol.nii.gz ' sprintf('TEMP/vol_%04dAffine.txt',j)]);
            temp_vol = load_untouch_nii(sprintf('ELASTIX_TEMP/vol_%04d/result.0.nii',j));
            ec_data.img(:,:,:,j) = temp_vol.img;
        end
        
        if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
            % READ ROTATION MATRIX
            f = fopen(sprintf('ELASTIX_TEMP/vol_%04d/TransformParameters.0.txt',j),'rt');
            while(isempty(strfind(fgetl(f),'CenterOfRotationPoint')))
                if(feof(f))
                    disp('Error reading transformation file');
                    fclose(f);
                    return
                end
            end
            useful_line = fgetl(f); % transform matrix
            useful_line = strrep(useful_line,'(','');
            useful_line = strrep(useful_line,')','');
            fclose(f);
            splitted_line = strsplit(useful_line,' ');
            RMatrix = [str2num(splitted_line{2}) str2num(splitted_line{3}) str2num(splitted_line{4})
                str2num(splitted_line{5}) str2num(splitted_line{6}) str2num(splitted_line{7})
                str2num(splitted_line{8}) str2num(splitted_line{9}) str2num(splitted_line{10})];
            ec_data.bvecs(j,:) = RMatrix*(dw_data.bvecs(j,:)');
            ec_data.bvecs(j,:) = ec_data.bvecs(j,:)/norm(ec_data.bvecs(j,:),2);
            yaw(j) = atan2d(RMatrix(2,1),RMatrix(1,1));
            pitch(j) = atan2d(-RMatrix(3,1),sqrt(RMatrix(3,2).^2+RMatrix(3,3).^2));
            roll(j) = atan2d(RMatrix(3,2),RMatrix(3,3));
            traslx(j) = str2num(splitted_line{11});
            trasly(j) = str2num(splitted_line{12});
            traslz(j) = str2num(splitted_line{13});
            determinant(j) = det(RMatrix);
        end
        
        if(exist('show_plot','var') > 0 && show_plot == 1)
            subplot(211)
            hold off
            plot(yaw,'r');
            hold on
            plot(pitch,'g');
            plot(roll,'b');
            title('Yaw pitch roll');
            ylabel('deg');
            legend({'RotX','RotY','RotZ'});
            subplot(212)
            hold off
            plot(traslx,'r');
            hold on
            plot(trasly,'g');
            plot(traslz,'b');
            ylabel('mm');
            title('Traslations X-Y-Z');
            legend({'X','Y','Z'});
            drawnow
        end
    end
    save_untouch_nii(ec_data,[dataname '_elastix_ec.nii.gz']);
    if(exist('rotate_bvecs', 'var') > 0 && rotate_bvecs == 1)
        bvecs = ec_data.bvecs';
        save([dataname '_elastix_ec.bvec'],'bvecs', '-ascii');
    end
    
    if(exist('show_plot','var') > 0 && show_plot == 1)
        savefig([dataname '_elastix.fig']);
    end
    save([dataname '_elastix'],'yaw','roll','pitch','determinant','traslx','trasly','traslz');
    rmdir('ELASTIX_TEMP','s');
end
end