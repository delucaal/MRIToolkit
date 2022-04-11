%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



classdef DistProbSHTracker < SHTracker
    % This is a probabilistic tracking class based on SHTracker. It first
    % reconstructs the FOD peaks. At each step, a uniform probabablity is assigned
    % to each peak and weighted by its relative amplitude. The regridding is performed
    % on the FODs.
    
    properties (Access = public)
        number_of_iterations=1
        sample_fod_peaks=0
        weight_mode=1
        stepSize_sd;
        maxAngle_sd;
    end
    
    properties (Access = protected)
        peak_dirs, peak_vals
        peak_dirs_points, peak_vals_points
        s2i
        stepSize_mu;
        maxAngle_mu;
    end
    
    methods (Access = public)
        function this = DistProbSHTracker(v2w, op, nof_iters)
            if(nargin < 2)
                op = [];
            end
            this = this@SHTracker(v2w,op);
            if(nargin > 2)
                this.number_of_iterations = nof_iters;
            end
            disp('Using class DistProbSHTracker for tracking (peak-based probabilistic tractography)');
            rng('default')
            rng(123456)
        end
        
        function this = setNumberOfIterations(this,iters)
            this.number_of_iterations = iters;
            disp(['Running ' num2str(iters) ' probabilistic tractography iterations']);
        end
        
        function this = setParametersSD(this,stepSize_sd,angle_sd)
            if(this.stepSize == 0)
                warning('Please call setParametersSD only after calling setParameters');
            end
            this.stepSize_sd = stepSize_sd;
            this.maxAngle_sd = angle_sd;
        end
        
        function this = setParameters(this, stepSize, threshold, maxAngle, lengthRange)
            if stepSize > 0
                this.stepSize = single(stepSize);
                this.stepSize_mu = single(stepSize);
            else
                error('stepSize must be > 0');
            end
            if threshold >= 0
                this.threshold = single(threshold);
            else
                error('threshold must be >= 0');
            end
            if maxAngle > 0 && maxAngle <= 90
                this.maxAngle = single(maxAngle);
                this.maxAngle_mu = single(maxAngle);
            else
                error('maxAngle must be in [0 90]');
            end
            if lengthRange(1) >= stepSize
                if (lengthRange(1) <= lengthRange(2))
                    this.lengthRange = single(lengthRange);
                else
                    error('lengthRange(2) must be >= lengthRange(1)');
                end
            else
                error('lengthRange(1) must be >= stepSize');
            end
        end
                
        function [totalTract, totalTractVal] = track(this, opoint)
            totalTract = {};
            totalTractVal = {};
            
            opoint = single(opoint);
            if (isempty(this.stepSize) || isempty(this.threshold) || isempty(this.maxAngle) || isempty (this.lengthRange))
                error('set tracking parameters first, using setParameters');
            end
            if isempty(this.f)
                error('set input data first, using setData');
            end
            
            this.peak_dirs = cell(size(this.f(:,:,:,1)));
            this.peak_vals = cell(size(this.f(:,:,:,1)));
            
            for track_id=1:this.number_of_iterations
                
                this.update_tracking_constraints();
                
                % interpolate
                this.interpolate(opoint);
                
                % mask out NaN values
                mask = ~isnan(this.fun(1,:));
                opoint = opoint(:,mask);
                this.fun = this.fun(:,mask);
                if(size(this.fun,2) ~= size(this.s2i,2))
                    this.s2i = this.s2i(mask);
                    this.update_peaks();
                end
                
                % process function
                this.process();
                
                % determine initial track direction(s)
                [point, dir, val] = this.getInitDir(opoint);
                
                % mask out all small peaks
                mask = val > this.threshold;
                point = point(:,mask);
                dir = dir(:,mask);
                
                % repeat data for probabilistic tractography
                point = repmat(point,[1 this.iterations]);
                dir = repmat(dir,[1 this.iterations]);
                val = repmat(val,[1 this.iterations]);
                
                % Track in both directions
                if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in first direction'); end;
                [tract1, tractVal1] = this.trackOneDir(point, dir);
                if this.pb; progressbar('ready'); end;
                
                if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in second direction'); end;
                [tract2, tractVal2] = this.trackOneDir(point, -dir);
                if this.pb; progressbar('ready'); end;
                
                % join tracts from both directions
                tract2 = cellfun(@flipud,tract2,'UniformOutput',false);
                tractVal2 = cellfun(@flipud,tractVal2,'UniformOutput',false);
                tract = cell(size(tract2));
                tractVal = cell(size(tractVal2));
                for j = 1:size(tract2,2)
                    if ~isempty(tract2{j})
                        tract{j} = [tract2{j}; point(:,j)'];
                        tractVal{j} = [tractVal2{j}; val(:,j)'];
                        if ~isempty(tract1{j})
                            tract{j} = [tract{j}; tract1{j}];
                            tractVal{j} = [tractVal{j}; tractVal1{j}];
                        end
                    else
                        if ~isempty(tract1{j})
                            tract{j} = [point(:,j)'; tract1{j}];
                            tractVal{j} = [val(:,j)'; tractVal1{j}];
                        end
                    end
                end
                
                % enforce length limitations
                %             maska = cellfun('size',tract,1)*this.stepSize >= this.lengthRange(1);
                %             maskb = cellfun('size',tract,1)*this.stepSize <= this.lengthRange(2);
                %             mask = maska & maskb;
                %             tract = tract(mask);
                %             tractVal = tractVal(mask);
                mask = cellfun('size',tract,1)>1;
                tract = tract(mask);
                tractVal = tractVal(mask);
                
                totalTract = [totalTract tract];
                totalTractVal = [totalTractVal tractVal];
            end
        end
        
    end
    
    methods (Access = protected)
        function process(this)
            if ~isempty(this.op);
                this.fun = this.op.process(this.fun);
            end
        end
        
        function [dir, val, angle] = getDir(this, prevDir)
            if(this.sample_fod_peaks == 1)
                [dir,val] = this.getSampledDir(prevDir);
            else
                [dir,val] = getDir@SHTracker(this,prevDir);
            end
            angle = (180/pi)*real(acos(abs(sum(prevDir.*dir,1))));
        end
        
        function [dir,val] = getSampledDir(this,prevDir)
            %             [adirs,avals] = SHPrecomp.all_peaks(this.fun,this.threshold,4);
            adirs = this.peak_dirs_points;
            avals = this.peak_vals_points;
            
            dir = zeros(3,length(avals));
            val = zeros(1,length(avals));
            
            for ij=1:length(avals)
                avals_t = avals{ij};
                if(isempty(avals_t))
                    continue
                end
                for ak=1:length(avals_t)
                    ddir = adirs{ij}(:,ak);
                    ang =  (180/pi)*real(acos(abs(sum(ddir.*prevDir(:,ij),1))));
                    if(ang > this.maxAngle)
                        avals_t(ak) = 0;
                    end
                end
                
                if(this.weight_mode == 0)
                    % Probability weighted by FOD amplitude
                    chances = avals_t/sum(avals_t);
                elseif(this.weight_mode == 1)
                    % Equal probability for all peaks
                    chances = ones(size(avals_t));
                    chances(avals_t == 0) = 0;
                    chances = chances/sum(chances);
                else
                    error('Unsupported weight mode');
                end
                chances_slot = randi(1000,1,length(avals_t)) .* chances;
                [~,IX] = max(chances_slot);
                val(ij) = avals_t(IX);
                dir(:,ij) = adirs{ij}(:,IX);
            end
            %             toc
        end
        
        function update_peaks(this)
            this.peak_dirs_points = this.peak_dirs(this.s2i);
            this.peak_vals_points = this.peak_vals(this.s2i);
            empty = find(cellfun('size', this.peak_vals_points, 2)  == 0);
            subset = this.s2i(empty);
            [this.peak_dirs_points(empty),this.peak_vals_points(empty)]  = SHPrecomp.all_peaks(this.fun(:,empty),this.threshold,4);
            this.peak_dirs(subset) = this.peak_dirs_points(empty);
            this.peak_vals(subset) = this.peak_vals_points(empty);
        end
        
        function interpolate(this, point)
            point(4,:) = 1; voxel = this.v2w\point; voxel = voxel(1:3,:);
            this.fun = mex_interp(this.f, voxel);
            voxel = round(voxel);
            voxel(voxel == 0) = 1;
            voxel(1,voxel(1,:) > size(this.f,1)) = size(this.f,1);
            voxel(2,voxel(2,:) > size(this.f,2)) = size(this.f,2);
            voxel(3,voxel(3,:) > size(this.f,3)) = size(this.f,3);
            this.s2i = sub2ind(size(this.f(:,:,:,1)),voxel(1,:),voxel(2,:),voxel(3,:));
            if(size(this.s2i,2) ~= size(this.fun,2))
                warning('this will crash');
            end
            this.update_peaks();
        end
        
        function mask = mask_nan(this)
            mask = ~isnan(this.fun(1,:));
            this.fun = this.fun(:,mask);
            this.s2i = this.s2i(mask);
            this.peak_dirs_points = this.peak_dirs_points(mask);
            this.peak_vals_points = this.peak_vals_points(mask);
        end
        
        function this = update_tracking_constraints(this)
            if(this.stepSize_sd ~= 0)
                this.stepSize = eps+max(normrnd(this.stepSize_mu,this.stepSize_sd),0);
            end
            if(this.maxAngle_sd ~= 0)
                this.maxAngle = eps+max(normrnd(this.maxAngle_mu,this.maxAngle_sd),0);
            end
        end
        
    end
end
