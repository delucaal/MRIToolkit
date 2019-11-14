%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function varargout = DW_Elastix_Register_GUI(varargin)
% DW_ELASTIX_REGISTER_GUI MATLAB code for DW_Elastix_Register_GUI.fig
%      DW_ELASTIX_REGISTER_GUI, by itself, creates a new DW_ELASTIX_REGISTER_GUI or raises the existing
%      singleton*.
%
%      H = DW_ELASTIX_REGISTER_GUI returns the handle to a new DW_ELASTIX_REGISTER_GUI or the handle to
%      the existing singleton*.
%
%      DW_ELASTIX_REGISTER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DW_ELASTIX_REGISTER_GUI.M with the given input arguments.
%
%      DW_ELASTIX_REGISTER_GUI('Property','Value',...) creates a new DW_ELASTIX_REGISTER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DW_Elastix_Register_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DW_Elastix_Register_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DW_Elastix_Register_GUI

% Last Modified by GUIDE v2.5 28-Mar-2018 09:47:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DW_Elastix_Register_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DW_Elastix_Register_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DW_Elastix_Register_GUI is made visible.
function DW_Elastix_Register_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DW_Elastix_Register_GUI (see VARARGIN)

% Choose default command line output for DW_Elastix_Register_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DW_Elastix_Register_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DW_Elastix_Register_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in availableListBox.
function availableListBox_Callback(hObject, eventdata, handles)
% hObject    handle to availableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns availableListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from availableListBox


% --- Executes during object creation, after setting all properties.
function availableListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to availableListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% RETRIEVE PATH TO REGISTRATION FILES
example = which('elastix_affine_klein_linterp.txt');
if(strfind(example,'/') > 0)
    path = strsplit(example,'/');
else
    path = strsplit(example,'\');
end
fpath = [];
for j=1:length(path)-1
    if(isempty(path{j}))
        path{j} = '/';
    end
    fpath = fullfile(fpath,path{j});
end
files = dir(fullfile(fpath,'*.txt'));
regNames = cell(1,length(files));
regFiles = cell(1,length(files));
for j=1:length(regNames)
   regNames{j} = files(j).name; 
   regFiles{j} = fullfile(fpath,files(j).name);
end
handles.regNames = regNames;
handles.regFiles = regFiles;
handles.AddedIndexes = [];
set(hObject,'String',regNames);
guidata(hObject,handles);


% --- Executes on button press in addTransformationButton.
function addTransformationButton_Callback(hObject, eventdata, handles)
% hObject    handle to addTransformationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
abox = findobj('Tag','availableListBox');
rbox = findobj('Tag','registrationListBox');
cidx = get(abox,'Value');
handles.AddedIndexes = union(handles.AddedIndexes,cidx,'stable');
set(rbox,'String',handles.regNames(handles.AddedIndexes));
set(rbox,'Value',1);
guidata(hObject,handles);


% --- Executes on button press in removeTransformationButton.
function removeTransformationButton_Callback(hObject, eventdata, handles)
% hObject    handle to removeTransformationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rbox = findobj('Tag','registrationListBox');
cidx = get(rbox,'Value');
set(rbox,'Value',1);
handles.AddedIndexes = setdiff(handles.AddedIndexes,handles.AddedIndexes(cidx));
set(rbox,'String',handles.regNames(handles.AddedIndexes));
guidata(hObject,handles);


% --- Executes on selection change in registrationListBox.
function registrationListBox_Callback(hObject, eventdata, handles)
% hObject    handle to registrationListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns registrationListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from registrationListBox


% --- Executes during object creation, after setting all properties.
function registrationListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to registrationListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function movingBox_Callback(hObject, eventdata, handles)
% hObject    handle to movingBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movingBox as text
%        str2double(get(hObject,'String')) returns contents of movingBox as a double


% --- Executes during object creation, after setting all properties.
function movingBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movingBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in movingChooseButton.
function movingChooseButton_Callback(hObject, eventdata, handles)
% hObject    handle to movingChooseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uigetfile('*.nii;*.nii.gz;*.hdr','Nifti files');
set(findobj('Tag','movingBox'),'String',[fp fn]);
handles.movingFile = [fp fn];
guidata(hObject,handles);



function fixedBox_Callback(hObject, eventdata, handles)
% hObject    handle to fixedBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedBox as text
%        str2double(get(hObject,'String')) returns contents of fixedBox as a double


% --- Executes during object creation, after setting all properties.
function fixedBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fixedChooseButton.
function fixedChooseButton_Callback(hObject, eventdata, handles)
% hObject    handle to fixedChooseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uigetfile('*.nii;*.nii.gz;*.hdr','Nifti files');
set(findobj('Tag','fixedBox'),'String',[fp fn]);
handles.fixedFile = [fp fn];
guidata(hObject,handles);


function outputDirBox_Callback(hObject, eventdata, handles)
% hObject    handle to outputDirBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputDirBox as text
%        str2double(get(hObject,'String')) returns contents of outputDirBox as a double


% --- Executes during object creation, after setting all properties.
function outputDirBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputDirBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputChooseButton.
function outputChooseButton_Callback(hObject, eventdata, handles)
% hObject    handle to outputChooseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fp = uigetdir('Output Directory');
set(findobj('Tag','outputDirBox'),'String',fp);
handles.outputDir = fp;
guidata(hObject,handles);


function outputFilenameBox_Callback(hObject, eventdata, handles)
% hObject    handle to outputFilenameBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputFilenameBox as text
%        str2double(get(hObject,'String')) returns contents of outputFilenameBox as a double


% --- Executes during object creation, after setting all properties.
function outputFilenameBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputFilenameBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputFilenameChooseButton.
function outputFilenameChooseButton_Callback(hObject, eventdata, handles)
% hObject    handle to outputFilenameChooseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uiputfile('*.gz;','Nifti files');
set(findobj('Tag','outputFilenameBox'),'String',[fp fn]);
handles.outputFilename = [fp fn];
guidata(hObject,handles);


% --- Executes on button press in registerButton.
function registerButton_Callback(hObject, eventdata, handles)
% hObject    handle to registerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
use_mask_checkbox = findobj('Tag','useMaskCheckbox');
use_mask = get(use_mask_checkbox,'Value');
if(use_mask == 0)
    fMask = '';
    mMask = '';
else
    fMask = handles.fixedMaskFile;
    mMask = handles.movingMaskFile;
end
DW_Elastix_Register(handles.movingFile,handles.fixedFile,handles.regFiles(handles.AddedIndexes),handles.outputDir,mMask,fMask,handles.outputFilename,handles);



function movingMaskBox_Callback(hObject, eventdata, handles)
% hObject    handle to movingMaskBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movingMaskBox as text
%        str2double(get(hObject,'String')) returns contents of movingMaskBox as a double


% --- Executes during object creation, after setting all properties.
function movingMaskBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movingMaskBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in movingMaskButton.
function movingMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to movingMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uigetfile('*.nii;*.nii.gz;*.hdr','Nifti files');
set(findobj('Tag','movingMaskBox'),'String',[fp fn]);
handles.movingMaskFile = [fp fn];
guidata(hObject,handles);


function fixedMaskBox_Callback(hObject, eventdata, handles)
% hObject    handle to fixedMaskBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of fixedMaskBox as text
%        str2double(get(hObject,'String')) returns contents of fixedMaskBox as a double


% --- Executes during object creation, after setting all properties.
function fixedMaskBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedMaskBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fixedMaskButton.
function fixedMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to fixedMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uigetfile('*.nii;*.nii.gz;*.hdr','Nifti files');
set(findobj('Tag','fixedMaskBox'),'String',[fp fn]);
handles.fixedMaskFile = [fp fn];
guidata(hObject,handles);

% --- Executes on button press in useMaskCheckbox.
function useMaskCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to useMaskCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useMaskCheckbox

function DW_Elastix_Register(moving,fixed,regparlist,output_dir,mMask,fMask,outfilename,handles)
%     elastix_path = '/Users/alb/Desktop/M_Code_ExploreDTI_v4.8.6/Source/MD_cor_E/macOSX64/';
%     elastix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'elastix_Mac_64'];
    elastix_cmd = handles.elastixLocation;
    sentence = [elastix_cmd ' -m ' moving ' -f ' fixed ' -out ' output_dir];
    if(exist('mMask','var') > 0 && ~isempty(mMask))
       sentence = [sentence ' -mMask ' mMask];
    end
    if(exist('fMask','var') > 0 && ~isempty(fMask))
       sentence = [sentence ' -fMask ' fMask];
    end
    for j=1:length(regparlist)
       sentence = [sentence ' -p ' regparlist{j}]; 
    end
    if(~exist(output_dir,'dir'))
        mkdir(output_dir)
    end
    system(sentence);
    if(exist('outfilename','var') > 0)
        result = dir([output_dir '/result*.nii.gz']);
        if(isempty(result))
           tp = dir([output_dir '/TransformParameters.*.txt']);
           DW_Elastix_Transform(moving,outfilename,[output_dir '/' tp(end).name],handles); 
        else
           if(length(result) == 1)  
               try
                  copyfile(fullfile(output_dir,result.name),outfilename);
               catch
                   system(['cp ' fullfile(output_dir,result.name) ' ' outfilename]);
               end
           else
               copyfile(fullfile(output_dir,result(end).name),outfilename);
               system(['cp ' fullfile(output_dir,result(end).name) ' ' outfilename]);
           end
        end
    end


function elastixLocation_Callback(hObject, eventdata, handles)
% hObject    handle to elastixLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elastixLocation as text
%        str2double(get(hObject,'String')) returns contents of elastixLocation as a double


% --- Executes during object creation, after setting all properties.
function elastixLocation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elastixLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [fn, fp] = uigetfile('*.exe','Elastix executable');
fp = uigetdir(pwd,'Elastix folder');
elastix = dir(fullfile(fp,'*elastix*.exe'));
transformix = dir(fullfile(fp,'*transformix*.exe'));
set(findobj('Tag','elastixLocation'),'String',fp);
handles.elastixLocation = fullfile(fp,elastix.name);
handles.transformixLocation = fullfile(fp,transformix.name);
guidata(hObject,handles);

function DW_Elastix_Transform(filein,fileout,transform_parameters,handles)
    mkdir('temp_transformix');
    transformix_cmd = handles.transformixLocation;
%     elastix_path = '/Users/alb/Desktop/M_Code_ExploreDTI_v4.8.6/Source/MD_cor_E/macOSX64/';
%     transformix_cmd = ['LD_LIBRARY_PATH=' elastix_path ' ' elastix_path 'transformix_Mac_64'];
    header = load_untouch_header_only(filein);
    if(header.dime.dim(1) == 3)
        disp('3D file');
        sentence = [transformix_cmd ' -in ' filein ' -out temp_transformix/ -tp ' transform_parameters];
        system(sentence);
        try
            gzip('temp_transformix/*.nii');
        catch
        end
        copyfile('temp_transformix/result.nii.gz',fileout);
    elseif(header.dime.dim(1) == 4)
        disp('4D file');
        data = load_untouch_nii(filein);
        temp_data = data;
        temp_data.hdr.dime.dim(1) = 3;
        temp_data.hdr.dime.dim(5) = 1;
        system('mkdir temp_transformix_4d');
        for vol_id=1:header.dime.dim(5)
            temp_data.img = data.img(:,:,:,vol_id);
            save_untouch_nii(temp_data,'./temp_transformix_4d/tmp_elastix_transform.nii.gz');
            DW_Elastix_Transform('./temp_transformix_4d/tmp_elastix_transform.nii.gz',sprintf('./temp_transformix_4d/reg_tmp_reg_elastix_%.3d.nii.gz',vol_id),transform_parameters);
        end
        system(['fslmerge -t ' fileout ' ./temp_transformix_4d/reg_tmp_reg_elastix*.gz']);
        system('rm -rf temp_transformix_4d');
    end
    try
        delete('temp_transformix/*');
        rmdir('temp_transformix');
    catch
    end