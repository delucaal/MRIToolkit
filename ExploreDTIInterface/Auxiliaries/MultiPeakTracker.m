%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

classdef MultiPeakTracker < SHTracker
    properties (Access = protected)
        nlevels = 1;
    end
    
    methods (Access = protected)
       
        function fun = interpolate_local(this, point)
            point(4,:) = 1; voxel = this.v2w\point; voxel = voxel(1:3,:);
            fun = mex_interp(this.f, voxel);
        end
        
        function [outTract,outTractVal] = postprocess(this,tract,tractVal)
            disp('Running the multipeak postprocess');
            outTract = tract;
            outTractVal = tractVal;
            for iter = 1:this.nlevels
                disp(['Iter ' num2str(iter)]);
                for tract_id=1:length(outTract)
                    TP = outTract{tract_id};
                    VAL = outTractVal{tract_id};
                    fods = this.interpolate_local(TP');
                                        
                    [~,amps] = SHPrecomp.all_peaks(fods,this.threshold,4);
                    ndirs = cellfun('length',amps);
                    NM = find(ndirs > 1);
                    new_seeds = TP(NM,:)';

                    try
                        % interpolate
                        this.interpolate(new_seeds);

                        % mask out NaN values
                        mask = this.mask_nan();
                        new_seeds = new_seeds(:,mask);

                        % process function
                        this.process();

                        % determine initial track direction(s)
                        [point, dir, val] = this.getInitDir(new_seeds);

                        % mask out all small peaks
                        mask = val > this.threshold;
                        point = point(:,mask);
                        dir = dir(:,mask);
                        
                        [newTract, newTractVal] = this.trackOneDir(point,dir);
                        
                        duplicates = false(length(newTract),1);
                        for new_tract_id = 1:length(newTract)
                            TP2 = newTract{new_tract_id};
                            if(size(TP2,1) > 1)
                                P2C = TP2(2,:);
                                % Check if we went back from the seed
                                % points (using the point before the split)
                                for point_id=1:length(NM)
                                    if(NM(point_id) > 1)
                                        if(norm(P2C-TP(NM(point_id)-1,:),2) < 1.1*this.stepSize)
                                            duplicates(point_id) = true;
                                        end
                                    end
                                end
                            else
                                duplicates(new_tract_id) = true;
                            end                                
                        end
                        newTract(duplicates) = [];
                        newTractVal(duplicates) = [];                     
                        
                        for new_tract_id = 1:length(newTract)
                            NT = newTract{new_tract_id};
                            NTV = newTractVal{new_tract_id};
                            for split_point = 1:size(TP,1)
                               D1 = norm(TP(split_point,:)-NT(1,:));
                               if(D1 < 1.1*this.stepSize)
                                   break
                               end
                            end
                            % Check if a flip is needed
                            V1 = (TP(split_point,:)-NT(1,:));
                            V2 = (TP(split_point,:)-NT(2,:));
                            NV1 = norm(V1);
                            if(NV1 > 1.1*this.stepSize)
                                disp('Strange');
                            end
                            
                            if(norm(NT(1,:)-TP(split_point,:)) < ...
                                    norm(NT(end,:)-TP(split_point,:)))
                                outTract(end+1) = {[TP(1:split_point,:);NT]};
                                outTractVal(end+1) = {[VAL(1:split_point);NTV]}; 
                                outTract(end+1) = {[flip(TP(split_point:end,:),1);NT]};
                                outTractVal(end+1) = {[flip(VAL(split_point:end));NTV]}; 
                            else
                                outTract(end+1) = {[TP(1:split_point,:);flip(NT)]};
                                outTractVal(end+1) = {[VAL(1:split_point);flip(NTV)]}; 
                                outTract(end+1) = {[flip(TP(split_point,:),1);flip(NT)]};
                                outTractVal(end+1) = {[flip(VAL(split_point:end));flip(NTV)]}; 
                            end                            

                        end
                    catch
                        disp('Check');
                    end
                end
            end
            
            % Finally, apply the length constraints
            lengths = cellfun('size',outTract,1)*this.stepSize;
            good_tracts = lengths >= this.lengthRange(1) & lengths <= this.lengthRange(2);
            outTract = outTract(good_tracts);
            outTractVal = outTractVal(good_tracts);
            
        end
        
    end
    
end
