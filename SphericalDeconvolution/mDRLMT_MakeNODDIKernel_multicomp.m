function [bmat,LRKernel,HRKernel] = mDRLMT_MakeNODDIKernel_multicomp(data,nreconstruction_vertices,lr_nreconstruction_vertices,noddi_values,isoDs,shell_data)
addpath(genpath('NODDI_toolbox_v1.0.1'));

shells = single(unique(int16(round(data.bvals)/1)*1)); % automatic detection of number of shells.
% This shell splitting works only for b-value spaced more than 100. To be fixed for other datasets.
ndirections = zeros(length(shells),1);
for ij=1:length(shells)
    ndirections(ij) = sum(abs(data.bvals-shells(ij))<1);
end

bvals = data.bvals';
save('temp_bvals.bval','bvals','-ascii');
bvecs = data.bvecs';
% bvecs(:,1) = -bvecs(:,1);
% bvecs(:,3) = -bvecs(:,3);
save('temp_bvecs.bvec','bvecs','-ascii');
clear bvals bvecs
NODDI_protocol = FSL2Protocol('temp_bvals.bval','temp_bvecs.bvec');
delete('temp_bvals.bval');
delete('temp_bvecs.bvec');


if(shell_data > 0)
    bmat = cell(length(shells),1); % bmat = [bvecs bvals]
    Kernel = cell(length(shells),1);
    Kernel_LR = cell(length(shells),1);
else
    bmat{1} = [data.bvecs data.bvals];
    Kernel = cell(1,1);
    Kernel_LR = cell(1,1);
end

% nreconstruction_vertices = 300; % 20/12/2017
nlr_vert = 0;
nhr_vert = sum(nreconstruction_vertices);

for ij=1:length(nreconstruction_vertices)
    super_scheme{ij} = gen_scheme(nreconstruction_vertices(ij),4); % the reconstruction scheme. Change 300 to any number
    [phi{ij}, theta{ij}] = cart2sph(super_scheme{ij}.vert(:,1),super_scheme{ij}.vert(:,2),super_scheme{ij}.vert(:,3)); % polar decomposition
    lr_scheme{ij} = gen_scheme(min(length(data.bvals),lr_nreconstruction_vertices(ij)),4);
    [phi_LR{ij}, theta_LR{ij}] = cart2sph(lr_scheme{ij}.vert(:,1),lr_scheme{ij}.vert(:,2),lr_scheme{ij}.vert(:,3));
    
    nlr_vert = nlr_vert + size(lr_scheme{ij}.vert,1);
end

HRKernel = zeros(sum(ndirections),nhr_vert+length(isoDs));
LRKernel = zeros(sum(ndirections),nlr_vert+length(isoDs));
% S = cell(length(shells),1);
% Ssim = zeros(sum(ndirections),1); % Just a simulated signal for internal testing

% noddi_values = [1 1.7E-9 3.5 0 3E-9 1];% x is the list of model parameters in SI units:
% % x(1) is the volume fraction of the intracellular space.
% % x(2) is the free diffusivity of the material inside and outside the cylinders.
% % x(3) is the concentration parameter of the Watson's distribution.
% % x(4) is the volume fraction of the isotropic compartment.
% % x(5) is the diffusivity of the isotropic compartment.
% % x(6) is the measurement at b=0.;

if(shell_data > 0)
    
    bvals = zeros(sum(ndirections),1);
    
    index = 1;
    
    for ij=1:length(shells)
        bmat{ij} = zeros(ndirections(ij),4);
        bmat{ij}(:,4) = shells(ij);
        bmat{ij}(:,1:3) = data.bvecs(index:index+ndirections(ij)-1,:);
        
        % Here the deconvolution dictionary is actually built
        
        % ANISOTROPIC PART
        hr_index_columns = 0;
        lr_index_columns = 0;
        
        Kernel{ij} = zeros(ndirections(ij),nhr_vert);
        Kernel_LR{ij} = zeros(ndirections(ij),nlr_vert);
        for aniso_comp = 1:length(nreconstruction_vertices)
            % HR
            for i=1:length(phi{aniso_comp})
                fibredir = super_scheme{aniso_comp}.vert(i,:)';
                E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
                Kernel{ij}(:,i+hr_index_columns) = E;
            end

            % LR
            for i=1:length(phi_LR{aniso_comp})
                fibredir = lr_scheme{aniso_comp}.vert(i,:)';
                E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
                Kernel_LR{ij}(:,i+lr_index_columns) = E;
            end
            
            hr_index_columns = hr_index_columns + size(super_scheme{aniso_comp}.vert,1);
            lr_index_columns = hr_index_columns + size(lr_scheme{aniso_comp}.vert,1);
        end
        
        % ISOTROPIC PART
        % HR
        for l=1:length(isoDs)
            Kernel{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end        
        
        % LR
        for l=1:length(isoDs)
            Kernel_LR{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end
        
        % Build a linear dictionary with the subparts
        bvals(index:index+ndirections(ij)-1) = bmat{ij}(:,4);
        HRKernel(index:index+ndirections(ij)-1,:) = Kernel{ij};
        LRKernel(index:index+ndirections(ij)-1,:) = Kernel_LR{ij};
        
        index = index+ndirections(ij);
    end
    
else
    N = length(data.bvals);
    bmat{1} = zeros(N,4);
    bmat{1}(:,1:3) = data.bvecs;
    bmat{1}(:,4) = data.bvals;

    % ANISOTROPIC PART
    hr_index_columns = 0;
    lr_index_columns = 0;
    for aniso_comp = 1:length(nreconstruction_vertices)
        
        % HR
        for i=1:length(phi{aniso_comp})
            fibredir = super_scheme{aniso_comp}.vert(i,:)';
            E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
            HRKernel(:,i+hr_index_columns) = E;
        end
        %LR
        for i=1:length(phi_LR{aniso_comp})
            fibredir = lr_scheme{aniso_comp}.vert(i,:)';
            E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
            LRKernel(:,i+lr_index_columns) = E;
        end
        
        hr_index_columns = hr_index_columns + size(super_scheme{aniso_comp}.vert,1);
        lr_index_columns = lr_index_columns + size(lr_scheme{aniso_comp}.vert,1);        
    end
    % ISOTROPIC PART
    % HR
    for l=1:length(isoDs)
        HRKernel(:,hr_index_columns+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
    % LR
    for l=1:length(isoDs)
        LRKernel(:,lr_index_columns+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end    
end

end