classdef MultiFieldTracker < handle
    %
    % Abstract Tractography class
    %
    % see also DTITracker, FODTracker
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (Access = protected)
        f;
        f2;
        v2w;
        stepSize;
        threshold;
        threshold2 = 0.3;
        maxAngle;
        maxAngle2 = 30;
        lengthRange;
        fun;
        fun2;
        iterations = 1;
        pb = false;
    end
    
    methods (Access = public)
        function this = MultiFieldTracker(v2w)
            this.v2w = single(v2w);
        end
        
        function this = setData(this, f, f2)
            this.f = single(f);
            this.f2 = single(f2);
        end

        function this = setMaxAngle2(this, angle2)
            this.maxAngle2 = angle2;
        end
        
        function this = setParameters(this, stepSize, threshold, maxAngle, lengthRange)
            if stepSize > 0
                this.stepSize = single(stepSize);
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
                if(isnan(this.maxAngle2))
                    this.maxAngle2 = this.maxAngle;
                end
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
        
        function this = setProgressbar(this, val)
            this.pb = val;
        end
        
        function [tract, tractVal] = track(this, point)
            point = single(point);
            if (isempty(this.stepSize) || isempty(this.threshold) || isempty(this.maxAngle) || isempty (this.lengthRange))
                error('set tracking parameters first, using setParameters');
            end
            if isempty(this.f)
                error('set input data first, using setData');
            end
            
            % interpolate
            this.interpolate(point);
            
            % mask out NaN values
            mask1 = ~isnan(this.fun(1,:));
            mask2 = ~isnan(this.fun2(1,:));
            point1 = point(:,mask1);
            point2 = point(:,mask2);
            this.fun = this.fun(:,mask1);
            this.fun2 = this.fun2(:,mask2);

            % process function
            this.process();
            
            % determine initial track direction(s)
            [point1, dir, val] = this.getInitDir(this.fun,point1);
            [point2, dir2, val2] = this.getInitDir(this.fun2,point2);            

            % mask out all small peaks
            mask1 = val > this.threshold;
            point1 = point1(:,mask1);
            dir = dir(:,mask1);

            % mask out all small peaks -2
            mask2 = val2 > this.threshold;
            point2 = point2(:,mask2);
            dir2 = dir2(:,mask2);            
            
            which_fod = cat(2,ones(1,size(point1,2)),2*ones(1,size(point2,2)));
            dir = cat(2,dir,dir2);
            val = cat(2,val,val2);
            point = cat(2,point1,point2);

            % repeat data for probabilistic tractography
            point = repmat(point,[1 this.iterations]);
            dir = repmat(dir,[1 this.iterations]);
            val = repmat(val,[1 this.iterations]);
            
            % Track in both directions
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in first direction'); end;
            [tract1, tractVal1] = this.trackOneDir(point, dir, which_fod);
            if this.pb; progressbar('ready'); end;
            
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in second direction'); end;
            [tract2, tractVal2] = this.trackOneDir(point, -dir, which_fod);
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
        end
    end
    
    methods (Access = private)
        function [tract, tractVal] = trackOneDir(this, point, dir, which_fod)
            tract = cell(1,size(point,2));
            tractVal = cell(1,size(point,2));
            flist = 1:size(point,2);
            
            for it = 1:(this.lengthRange(2)/this.stepSize)
                if this.pb; progressbar(size(point,2),int2str(size(point,2))); end;
                % advance streamline
                try

                    point = point + this.stepSize .* dir;
                catch
                    disp('Size point:')
                    disp(num2str(size(point)))
                    disp('Size step:')
                    disp(num2str(size(this.stepSize)))
                    disp('Size dir:')
                    disp(num2str(size(dir)))
                    continue
%                     error('The error...')
                end
                
                this.update_tracking_constraints();
                
                % interpolate
                this.interpolate(point);
                
                % mask out NaN values
                mask = ~isnan(this.fun(1,:)) | ~isnan(this.fun2(1,:));
                if(sum(mask) == 0)
                    disp('Will stop')
                end
                point = point(:,mask);
                dir = dir(:,mask);
                this.fun = this.fun(:,mask);
                this.fun2 = this.fun2(:,mask);
                flist = flist(mask);
                which_fod = which_fod(mask);

                % process function
                this.process();
                
                % get new direction
                [newDir, val, angle] = this.getDir(this.fun,dir);
                [newDir2, val2, angle2] = this.getDir(this.fun2,dir);

                % mask out small peaks
                mask = val > this.threshold | val2 > this.threshold2;
                if(sum(mask) == 0)
                    disp('Will stop')
                end
                point = point(:,mask);
                dir = dir(:,mask);                
                newDir = newDir(:,mask);
                newDir2 = newDir2(:,mask);
                flist = flist(mask);
                angle = angle(mask);
                angle2 = angle2(mask);
                val = val(:,mask);
                val2 = val2(:,mask);
                which_fod = which_fod(mask);

                % mask out large angles
                % Now check if you can jump into the second field
                mask1_angle = angle < this.maxAngle;
                mask1_amp = val > this.threshold;
                mask2_amp = val2 > this.threshold2; 
                mask2_angle = angle2 < this.maxAngle2;

                mask = true(size(mask1_angle));
                which_fod_new = NaN(size(which_fod));
                for ix=1:length(which_fod_new)
                    if(which_fod(ix) == 1)
                        % Currently following the first FOD
                        % Check if it can continue
                        if(mask1_amp(ix) > 0 && mask1_angle(ix) > 0)
                            % Yes, so let's leave it going
                            which_fod_new(ix) = 1;
                        else
                            % See if it can jump to the 2nd FOD
                            if(mask2_amp(ix) > 0)
                                which_fod_new(ix) = 2;
                                [d,v] = SHPrecomp.all_peaks(this.fun2(:,ix),this.threshold,4);
                                if(~isempty(v) && ~isempty(v{1}))
                                    newDir2(:,ix) = d{1}(:,1);
                                    val2(ix) = v{1}(1);
                                end
                            else
                                % Stop
                                which_fod_new(ix) = NaN;
                            end
                        end
                    elseif(which_fod(ix) == 2)
                        % Currently following the second FOD
                        % Check if it can continue
                        if(mask2_amp(ix) > 0 && mask2_angle(ix) > 0)
                            % Yes, so let's leave it going
                            which_fod_new(ix) = 2;
                        else
                            % See if it can jump to the 2nd FOD
                            if(mask1_amp(ix) > 0)
                                which_fod_new(ix) = 1;
                                [d,v] = SHPrecomp.all_peaks(this.fun(:,ix),this.threshold,4);
                                if(~isempty(v) && ~isempty(v{1}))
                                    newDir(:,ix) = d{1}(:,1);
                                    val(ix) = v{1}(1);
                                end
                            else
                                % Stop
                                which_fod_new(ix) = NaN;
                            end
                        end
                    end
                end
%                 if(sum(~isnan(angle2)) > 0)
%                     disp('Debug')
%                 end

                newDir(:,which_fod_new==2) = newDir2(:,which_fod_new==2);
                val(:,which_fod_new==2) = val2(:,which_fod_new==2);
                mask(which_fod_new==2) = 1;
                mask(isnan(which_fod_new)) = 0;
                which_fod = which_fod_new(mask);

                point = point(:,mask);
                dir = dir(:,mask);
                newDir = newDir(:,mask);
                flist = flist(mask);
                val = val(:,mask);
                
                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*newDir,1));
                
                % update dir
                dir = flipsign([1 1 1],:).*newDir;
                
                % stop if we are out of points
                if isempty(point)
                    disp('Out of points')
                    break
                end
                
                % add points to the tracts
                for i=1:length(flist)
                    tract{flist(i)}(it,:) = point(:,i);
                    tractVal{flist(i)}(it,:) = val(:,i);
                end
            end
        end
        
        function interpolate(this, point)
            point(4,:) = 1; voxel = this.v2w\point; voxel = voxel(1:3,:);
            this.fun = mex_interp(this.f, voxel);
            this.fun2 = mex_interp(this.f2, voxel);
        end
        
        function this = update_tracking_constraints(this)
            
        end
        
    end
    methods (Access = protected, Abstract = true)
        process(this);
        [point, dir, val, dir2, val2] = getInitDir(this, point);
        [dir, val, angle, dir2, val2, angle2] = getDir(this, prevDir);
    end
end