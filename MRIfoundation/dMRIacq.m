% MRIKaleidoscope MRI foundation class - A. De Luca
% Base diffusion MRI sequence class
% History:
% v1: 24/09/2018

classdef dMRIacq < MRIacq
    properties (Constant)
       ACQTYPE_SEDIFF = 'SEDIFF';
       ACQTYPE_GEDIFF = 'GEDIFF';
    end
    
    properties
       bvals % diffusion weightings (s/mm2)
       bvecs % diffusion gradients
       G % gradient strength of each volume (mT/m)
       delta % diffusion gradient width of each volume (s)
       DELTA % the time between gradients (s)
    end
    
    methods
        function obj = dMRIacq(obj_in)
           if(nargin > 0 && isa(obj_in,'MRIacq'))
               obj.TE = obj_in.TE;
               obj.TR = obj_in.TR;
               obj.TI = obj_in.TI;
               obj.AcqType = obj_in.AcqType;
               obj.AcqMode = obj_in.AcqMode;
               obj.hdr = obj_in.hdr;
               obj.img = obj_in.img;              
           end
           
           if(isfield(obj.hdr,'dime'))
            dim4 = obj.fourth_dimension();
           else
               dim4 = 1;
           end
            obj.TI = -1*ones(dim4,1);
            obj.TR = -1*ones(dim4,1);
            obj.TE = -1*ones(dim4,1);
            obj.AcqMode = -1*ones(dim4,1);
            obj.bvals = -1*ones(dim4,1);
            obj.bvecs = -1*ones(dim4,1);
            obj.G = -1*ones(dim4,1);
            obj.delta = -1*ones(dim4,1);
            obj.DELTA = -1*ones(dim4,1);
        end
        
        function obj = load_bvals_bvecs(obj,nii_file)
            fparts = strsplit(nii_file,'.');
            fname = fparts{1};
            for ij=2:length(fparts)
               if(strcmp(fparts{ij},'nii'))
                   break
               end
               fname = [fname '.' fparts{ij}];
            end
            obj.bvals = load([fname '.bval']);
            obj.bvals = obj.bvals';
            obj.bvecs = load([fname '.bval']);
            obj.bvecs = obj.bvecs';
        end
        
        % This function makes a handy structure of the main parameters, and
        % is used to ease the application of OptimizationMethods
        function prop_struct = properties_struct(obj)
           prop_struct = properties_struct@MRIacq(obj);
           prop_struct.bvals = obj.bvals;
           prop_struct.bvecs = obj.bvecs;
           prop_struct.G = obj.G;
           prop_struct.delta = obj.delta;
           prop_struct.DELTA = obj.DELTA;
        end
        
    end
    
    methods (Static)
        % Static data loader (without properties)
        % Input:
        % destination_string: the .nii file
        % Output:
        % obj: an MRIacq instance (with properties to be filled manually)
        function obj = load_volume(destination_string)
            obj = dMRIacq();
            tmp = load_untouch_nii(destination_string);
            obj.hdr = tmp.hdr;
            obj.img = single(tmp.img)*tmp.hdr.dime.scl_slope;
            
            dim4 = obj.fourth_dimension();
            obj.TI = -1*ones(dim4,1);
            obj.TR = -1*ones(dim4,1);
            obj.TE = -1*ones(dim4,1);
            obj.AcqMode = -1*ones(dim4,1);
            obj.bvals = -1*ones(dim4,1);
            obj.bvecs = -1*ones(dim4,1);
            obj.G = -1*ones(dim4,1);
            obj.delta = -1*ones(dim4,1);
            obj.DELTA = -1*ones(dim4,1);
        end
        % Static data loader (with properties)
        % Input:
        % destination_string: the .nii file
        % properties: a cell array with 2 string entries {'Prop','Value'}
        % Output:
        % obj: an MRIacq instance 
        function obj = load_volume_with_properties(destination_string,properties)
            obj = dMRIacq.load_volume(destination_string);
            nprops = length(properties)/2;
            for ij=1:nprops
               prop_name = properties{1+(ij-1)*2};
               if(~isprop(obj,prop_name))
                   error(['Not existing property ' prop_name]);
               end
               prop_value = properties{2+(ij-1)*2};
               if(length(prop_value) == 1 && isnumeric(prop_value))
                   prop_value = repmat(prop_value,obj.fourth_dimension(),1);
               end
               for vol_id=1:obj.fourth_dimension()
                   if(isnumeric(prop_value))
                       eval(['obj.' prop_name '(vol_id)=' num2str(prop_value(vol_id)) ';']);
                   else
                       eval(['obj.' prop_name '=''' prop_value ''';']);
                   end
               end
            end
        end
        
        % Static data loader (with properties)
        % Input:
        % destination_strings: a cell array of .nii files
        % properties: a cell array of cells with 2 string entries {'Prop','Value'}
        % Output:
        % obj: an MRIacq instance 
        function obj = load_volumes_with_properties(destination_strings,properties)
            obj = dMRIacq.load_volume_with_properties(destination_strings{1},properties{1});
            for ij=2:length(destination_strings)
               obj2 = dMRIacq.load_volume_with_properties(destination_strings{ij},properties{ij});
               obj = obj.stack_acquisitions(obj2);
            end
        end
        
    end
end