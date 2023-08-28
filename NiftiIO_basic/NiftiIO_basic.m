%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



classdef NiftiIO_basic < handle
   
    methods (Static)
        
        function ConvertDicomFolder2Nifti(dicom_folder,output_folder)
            global MRIToolkit;
            if(isfield(MRIToolkit,'dcm2niix') && exist(MRIToolkit.dcm2niix,'file'))
                cmd = [MRIToolkit.dcm2niix ' -o ' output_folder ' -b y -f %f_%p ' dicom_folder];
                status = system(cmd);
                if(status ~= 0)
                    warning('Something went wrong in the DICOM conversion');
                end
            else
                error('Cannot find dcm2niix or the location has not been specified in MRIToolkitDefineLocalVars.m');
            end
        end
        
        function OrganizeDataBIDSLike(dicom_folder,output_root)   
            global MRIToolkit;
            if(isempty(which('jsonread')))
                if(isfield(MRIToolkit,'JSONio') && exist(MRIToolkit.JSONio,'dir') > 0)
                    addpath(MRIToolkit.JSONio);
                else
                    error('This function needs the JSONio library to properly work. Configure MRIToolkitDefineLocalVars.m accordingly.');
                end
            end
            
            TempDirectory = fullfile(output_root,'TempProcessing');
            RawFolder = fullfile(output_root,'RawNii');
            OutputFolder = fullfile(output_root,'derivatives');

            if(exist(RawFolder,'dir') < 1)
                mkdir(RawFolder);
            end

            if(exist(OutputFolder,'dir') < 1)
                mkdir(OutputFolder);
            end

            if(exist(TempDirectory,'dir') < 1)
                mkdir(TempDirectory);
            end

            AcquisitionTypes = {};
            
            source_file = which('BIDSdefinitions.txt');
            if(isempty(source_file))
                disp('Cannot find BIDSdefinitions.txt in the path. Using built-in definitions');
                AcquisitionCategories = {
                    {'t1','anat'},...
                    {'t2','anat'},...
                    {'flair','anat'},...
                    {'dwi','dwi'},...
                    {'dti','dwi'},...
                    {'dmri','dwi'},...
                    {'fmri','func'},...
                    {'bold','func'}
                };
            else
                AcquisitionCategories = {};
                f = fopen(source_file,'rt');
                while(~feof(f))
                   l = fgetl(f);
                   if(contains(l,'#') || isempty(l))
                       continue
                   end
                   parts = strsplit(l,',');
                   AcquisitionCategories(end+1) = {parts};
                end
                fclose(f);
            end

            SubjectsFolder = dir(dicom_folder);
            for subj_id=1:length(SubjectsFolder)
                if(strcmp(SubjectsFolder(subj_id).name(1),'.'))
                    continue
                end
                delete(fullfile(TempDirectory,'*'));
                full_path = fullfile(dicom_folder,SubjectsFolder(subj_id).name);
                if(contains(SubjectsFolder(subj_id).name,'sub-') == false)
                    full_out_path = fullfile(RawFolder,['sub-' SubjectsFolder(subj_id).name]);
                else
                    full_out_path = fullfile(RawFolder,SubjectsFolder(subj_id).name);
                end
                npaths = length(dir([full_out_path '*']))+1;
                final_prefix = SubjectsFolder(subj_id).name;
                if(contains(SubjectsFolder(subj_id).name,'sub-') == false)
                    final_prefix = ['sub-' final_prefix];
                end
                if(contains(SubjectsFolder(subj_id).name,'-ses-') == false)
                    final_prefix = [final_prefix '-ses-' num2str(npaths)];
                    full_out_path = strcat(full_out_path,['-ses-' num2str(npaths)]);
                end
                
                NiftiIO_basic.ConvertDicomFolder2Nifti(full_path,TempDirectory);

                json_files = dir(fullfile(TempDirectory,'*.json'));
                if(~isempty(json_files))
                    for json_id=1:length(json_files)
                       file_content = jsonread(fullfile(TempDirectory,json_files(json_id).name)); 
                       seq_idx = StringInList(AcquisitionTypes,file_content.ProtocolName);
                       if(isempty(seq_idx))
                           AcquisitionTypes(end+1) = {file_content.ProtocolName};
                           seq_idx = length(AcquisitionTypes);
                       end
                       acq_category = AcquisitionTypes{seq_idx};
                       acq_sub_type = ForecastImgCategory(AcquisitionCategories,acq_category);
                       bvals_or_bvecs = 0;
                       if(strcmpi(acq_sub_type,'dwi'))
                           bvals_or_bvecs = 1;
                       else
                           gradients_info = dir(fullfile(TempDirectory,[json_files(json_id).name(1:end-5) '.bval']));               
                           if(~isempty(gradients_info))
                               bvals_or_bvecs = 1;
                               acq_sub_type = 'dwi';
                           end
                       end
                       if(exist(fullfile(full_out_path,acq_sub_type),'dir') < 1)
                           mkdir(fullfile(full_out_path,acq_sub_type));
                       end
%                        if(exist(fullfile(OutputFolder,acq_sub_type),'dir') < 1)
%                            mkdir(fullfile(OutputFolder,acq_sub_type));
%                        end
                       dest_prefix = fullfile(full_out_path,acq_sub_type,final_prefix);
                       dest_prefix = strcat(dest_prefix,['_' acq_category]);
                       existing_outputs = length(dir([dest_prefix '*.nii*']))+1;
                       dest_prefix = strcat(dest_prefix, ['_' num2str(existing_outputs)]);
                       copyfile(fullfile(TempDirectory,json_files(json_id).name),[dest_prefix '.json']);
                       corresponding_image = dir(fullfile(TempDirectory,[json_files(json_id).name(1:end-5) '.nii*']));
                       [~,~,extension] = fileparts(corresponding_image.name);
                       copyfile(fullfile(TempDirectory,corresponding_image.name),[dest_prefix extension]);
                       if(bvals_or_bvecs == 1)
                           bvals = fullfile(TempDirectory,[json_files(json_id).name(1:end-5) '.bval']);
                           bvecs = fullfile(TempDirectory,[json_files(json_id).name(1:end-5) '.bvec']);
                           try
                              copyfile(bvals,[dest_prefix '.bval']);
                              copyfile(bvecs,[dest_prefix '.bvec']);
                           catch
                           end
                       end
                    end
                end
            end
        end
 
        function description = ReadJSONDescription(varargin)
           % Read a .JSON file description
           % Input arguments:
           % json_file: the file to read
           
            json_file = GiveValueForName(varargin,'json_file');
            if(isempty(json_file))
                error('Missing mandatory argument json_file');
            end
           description = jsondecode(read_text_file(json_file));
        end
        
        function WriteJSONDescription(varargin)
            % Write a .JSON file describing the performed processing.
            % Input arguments:
            % props: a struct containing the properties to store to the
            % .json file. Can also be a file from which to copy most
            % parameters.
            % output: the output file excluding the .json extension
            global MRIToolkit
            
            if(isfield(MRIToolkit,'EnforceJSON'))
                if(MRIToolkit.EnforceJSON == false)
                    return
                end
            end
            
            test_struct.test = 'A';
            try
                jsonencode(test_struct);
            catch
                disp('There seems to be no JSON support. Disabling it.');
                MRIToolkit.EnforceJSON = false;
                return
            end
                
            jsonprops = GiveValueForName(varargin,'props');
            if(isempty(jsonprops))
                jsonprops = struct;
            elseif(ischar(jsonprops))
                jsonprops = jsondecode(read_text_file(jsonprops));
            end
            
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
                        
            CallerFunction = dbstack;
            jsonprops.ExecutedFunction = [CallerFunction(end).file ':' num2str(CallerFunction(end).line)];
            jsonprops.SoftwareVersion = MRIToolkit.version;
            
            jsonprops = jsonencode(jsonprops,'ConvertInfAndNaN',false);
            write_text_file([output '.json'],jsonprops); 
        end
        
    end
    
end

function acq_type = ForecastImgCategory(categories,img_name)
    acq_type = 'anat'; 
    for ij=1:length(categories)
        if(contains(img_name,categories{ij}{1},'IgnoreCase',true))
           acq_type = categories{ij}{2};
           return;
        end
    end
end

function idx = StringInList(the_list,the_string)
    idx = [];
    for ij=1:length(the_list)
        if(strcmpi(the_list{ij},the_string))
            idx = ij;
            return;
        end
    end
end

% Helper: finds a parameter by name when using varargin
function value = GiveValueForName(coptions,name)
value = [];
for ij=1:2:length(coptions)
    if(strcmpi(coptions{ij},name))
        value = coptions{ij+1};
        return
    end
end
end

function text = read_text_file(file)
    f = fopen(file,'rt');
    text = [];
    while(~feof(f))
        text = [text;fgetl(f)];
    end
    fclose(f);
end

function write_text_file(file,txt)
    f = fopen(file,'wt');
    txt = strrep(txt,'\n',newline);
    txt = strrep(txt,',"',[',' newline '"']);
    fprintf(f,'%s',txt);
    fclose(f);
end
