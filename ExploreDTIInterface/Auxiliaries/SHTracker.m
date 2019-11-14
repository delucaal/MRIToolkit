%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
classdef SHTracker < Tracker
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
%
    % SH Tractography class
    %
    % Usage:
    % tracker1 = SHTracker(v2w, sd/csd/qbi-object); % interpolation on dwi_sh followed by sd/csd/qbi-operation
    % tracker1.setParameters(stepSize, threshold, maxAngle, lengthRange);
    % tracker1.setData(dwi_sh);
    % tracker1.track(seedPoint);
    %
    % tracker2 = SHTracker(v2w); % interpolation directly on fod_sh
    % tracker2.setParameters(stepSize, threshold, maxAngle, lengthRange);
    % tracker2.setData(fod_sh);
    % tracker2.track(seedPoint);
    % 
    % see also Tracker, QBI, SD, CSD
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (Access = public)
        op;
    end
    
    methods (Access = public)
        function this = SHTracker(v2w, op)
            this = this@Tracker(v2w);
            if nargin > 1
                this.op = op;
                SHPrecomp.init(this.op.lmax);
            else
                this.op = [];
            end
        end
    
        function this = setData(this, f)
            this.f = single(f);
            if isempty(this.op)
                SHPrecomp.init(SH.n2lmax(size(f,4)));
            end
        end
    end
    
    methods (Access = protected)
        function process(this)
            if ~isempty(this.op);
                this.fun = this.op.process(this.fun);
            end
        end
        
        function [point, dir, val] = getInitDir(this, point)
            [dir, val] = SHPrecomp.all_peaks(this.fun, this.threshold, 4);
            c = cellfun('size', dir, 2);
            idx = [];
            for i = 1:size(c,2)
               idx = [idx ones(1,c(i))*i];
            end
            point = point(:,idx);
            
            dir = cellfun(@single,dir,'UniformOutput',false);
            val = cellfun(@single,val,'UniformOutput',false);
            
            dir = cell2mat(dir);
            val = cell2mat(val);
        end
        
        function [dir, val, angle] = getDir(this, prevDir)
            [dir, val] = SHPrecomp.peaks(this.fun, prevDir);
            angle = (180/pi)*real(acos(abs(sum(prevDir.*dir,1))));
        end
    end
end