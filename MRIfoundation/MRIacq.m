% MRIKaleidoscope MRI foundation class - A. De Luca
% Base MRI sequence class
% History:
% v1: 24/09/2018

classdef MRIacq
    properties
       TE % Echo time (s)
       TR % Repetition time (s)
       TI % Inversion time (s) (if any)
       AcqType % 'SE','GE','IR' (SpinEcho-GradientEcho-InversionRecovery)
       AcqMode % '2D-MS', 2 = '3D' (2D multi slice, 3D)
       hdr % Associated NIFTI header
       img % The ND volume data
    end
    
    properties (Constant)
       ACQTYPE_SE = 'SE'
       ACQTYPE_GE = 'GE'
       ACQTYPE_IR = 'IR'
       ACQMODE_2DMS = 1
       ACQMODE_3D = 0
    end
    
    methods
        % Initialize values
        function obj = MRIacq(varargin)
            obj.TE = -1;
            obj.TR = -1;
            obj.TI = -1;
            obj.AcqType = 'Invalid';
            obj.AcqMode = 'Invalid';
            obj.hdr.invalid = 1;
            obj.img = -1;
            for args_id=1:2:length(varargin)
                if(strcmp(varargin{args_id},'TE'))
                   obj.TE = varargin{args_id+1};
                elseif(strcmp(varargin{args_id},'TR'))
                   obj.TR = varargin{args_id+1};
                elseif(strcmp(varargin{args_id},'TI'))
                   obj.TI = varargin{args_id+1};
                elseif(strcmp(varargin{args_id},'img'))
                   obj.img = varargin{args_id+1};
                   obj.hdr.dime.dim = [ndims(obj.img) size(obj.img)];
                elseif(strcmp(varargin{args_id},'AcqType'))
                   obj.AcqType = varargin{args_id+1};
                elseif(strcmp(varargin{args_id},'AcqMode'))
                   obj.AcqMode = varargin{args_id+1};                   
                end
            end
        end
        
        % This function makes a handy structure of the main parameters, and
        % is used to ease the application of OptimizationMethods
        function prop_struct = properties_struct(obj)
           prop_struct.TE = obj.TE;
           prop_struct.TR = obj.TR;
           prop_struct.TI = obj.TI;
           prop_struct.AcqType = obj.AcqType;
           prop_struct.AcqMode = obj.AcqMode;
        end
        
        % General function to concatenate multiple acquisitions
        function obj = stack_acquisitions(obj1,obj2)
            obj = MRIacq();
            obj.TE = [obj1.TE;obj2.TE];
            obj.TR = [obj1.TR;obj2.TR];
            obj.TI = [obj1.TI;obj2.TI];
            if(~strcmp(obj1.AcqType,obj2.AcqType))
                error('Different acquisition types, unknown handling...');
            end
            obj.AcqType = obj1.AcqType;
            if(length(intersect(obj1.AcqMode,obj2.AcqMode)) ~= length(union(obj1.AcqMode,obj2.AcqMode)))
                warning('Different acquisition modes, continuing but expect errors.');
            end
            obj.AcqMode = [obj1.AcqMode;obj2.AcqMode];
            
            obj.hdr = obj1.hdr;
            [sx1,sy1,sz1,st1] = size(obj1.img);
            [sx2,sy2,sz2,st2] = size(obj2.img);
            if(sx1 ~= sx2 || sy1 ~= sy2 || sz1 ~= sz2)
                error('Different geometries, aborting.');
            end
            obj.img = zeros(sx1,sy1,sz1,st1+st2); 
            obj.img(:,:,:,1:st1) = obj1.img;
            obj.img(:,:,:,st1+1:st1+st2) = obj2.img;
            
            obj.hdr.dime.dim(5) = st1+st2;
            if(st1+st2 > 1)
                obj.hdr.dime.dim(1) = 4;
            else
                obj.hdr.dime.dim(1) = 3;
            end
        end
        
        % Prepare img for fast fit
        function obj = vectorize_data(obj)
            msize = obj.hdr.dime.dim(2:5);
            obj.img = reshape(obj.img,msize(1)*msize(2)*msize(3),msize(4));
        end
        
        % Recover original img dim
        function obj = reshape_data_to_original(obj)
            msize = obj.hdr.dime.dim(2:5);
            obj.img = reshape(obj.img,msize(1),msize(2),msize(3),msize(4));
        end
        
        % Discard specific volumes
        function obj = discard_volumes(obj,indices)
           obj.img(:,:,:,indices) = [];
           obj.hdr.dime.dim(5) = size(obj.img,4);
           obj.TI(indices) = [];
           obj.TE(indices) = [];
           obj.TR(indices) = [];
           if(~ismatrix(obj.AcqMode))
               obj.AcqMode(indices) = [];
           end
        end
        
        % Return 3D matrix size
        function mat = matrix_size(obj)
            mat = obj.hdr.dime.dim(2:4);
        end
        
        % Return only 4D matrix dimension
        function nd = fourth_dimension(obj)
            nd = obj.hdr.dime.dim(5);
        end
        
        % Sort img and properties according to property_name (must be numeric)
        % order: 'ascend' or 'descend'
        function obj = sort_by(obj1,property_name,order)
           eval(['my_prop = obj1.' property_name ';']);
           [~,IX] = sort(my_prop,order);
           obj = obj1;
           obj.TE = obj.TE(IX);
           obj.TR = obj.TR(IX);
           obj.TI = obj.TI(IX);
%            obj.AcqMode = obj.AcqMode(IX);
           obj.img = obj.img(:,:,:,IX);
        end
        
        % Normalize by 4th dimension max
        function obj = normalize_4d(obj)
           V = max(obj.img,[],4); 
           for ij=1:obj.fourth_dimension()
              obj.img(:,:,:,ij) = obj.img(:,:,:,ij)./(V+eps); 
           end
        end
        
    end
    
    methods (Static)
        % Static data loader (without properties)
        % Input:
        % destination_string: the .nii file
        % Output:
        % obj: an MRIacq instance (with properties to be filled manually)
        function obj = load_volume(destination_string)
            obj = MRIacq();
            tmp = load_untouch_nii(destination_string);
            obj.hdr = tmp.hdr;
            obj.img = single(tmp.img)*tmp.hdr.dime.scl_slope;
            
            dim4 = obj.fourth_dimension();
            obj.TI = -1*ones(dim4,1);
            obj.TR = -1*ones(dim4,1);
            obj.TE = -1*ones(dim4,1);
            obj.AcqMode = -1*ones(dim4,1);
        end
        
        % Static data loader (with properties)
        % Input:
        % destination_string: the .nii file
        % properties: a cell array with 2 string entries {'Prop','Value'}
        % Output:
        % obj: an MRIacq instance 
        function obj = load_volume_with_properties(destination_string,properties)
            obj = MRIacq.load_volume(destination_string);
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
            obj = MRIacq.load_volume_with_properties(destination_strings{1},properties{1});
            for ij=2:length(destination_strings)
               obj2 = MRIacq.load_volume_with_properties(destination_strings{ij},properties{ij});
               obj = obj.stack_acquisitions(obj2);
            end
        end
        
    end
end
