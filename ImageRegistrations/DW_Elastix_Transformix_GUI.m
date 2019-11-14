%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
function varargout = DW_Elastix_Transformix_GUI(varargin)
% DW_ELASTIX_TRANSFORMIX_GUI MATLAB code for DW_Elastix_Transformix_GUI.fig
%      DW_ELASTIX_TRANSFORMIX_GUI, by itself, creates a new DW_ELASTIX_TRANSFORMIX_GUI or raises the existing
%      singleton*.
%
%      H = DW_ELASTIX_TRANSFORMIX_GUI returns the handle to a new DW_ELASTIX_TRANSFORMIX_GUI or the handle to
%      the existing singleton*.
%
%      DW_ELASTIX_TRANSFORMIX_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DW_ELASTIX_TRANSFORMIX_GUI.M with the given input arguments.
%
%      DW_ELASTIX_TRANSFORMIX_GUI('Property','Value',...) creates a new DW_ELASTIX_TRANSFORMIX_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DW_Elastix_Transformix_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DW_Elastix_Transformix_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DW_Elastix_Transformix_GUI

% Last Modified by GUIDE v2.5 07-Dec-2015 16:45:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DW_Elastix_Transformix_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DW_Elastix_Transformix_GUI_OutputFcn, ...
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


% --- Executes just before DW_Elastix_Transformix_GUI is made visible.
function DW_Elastix_Transformix_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DW_Elastix_Transformix_GUI (see VARARGIN)

% Choose default command line output for DW_Elastix_Transformix_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DW_Elastix_Transformix_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DW_Elastix_Transformix_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


function transformParametersBox_Callback(hObject, eventdata, handles)
% hObject    handle to transformParametersBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transformParametersBox as text
%        str2double(get(hObject,'String')) returns contents of transformParametersBox as a double


% --- Executes during object creation, after setting all properties.
function transformParametersBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transformParametersBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in transformChooseButton.
function transformChooseButton_Callback(hObject, eventdata, handles)
% hObject    handle to transformChooseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, fp] = uigetfile('*.txt;','Transform Parameters');
set(findobj('Tag','transformParametersBox'),'String',[fp fn]);
handles.transformParameters = [fp fn];
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
[fn, fp] = uiputfile('*.gz','Nifti files');
set(findobj('Tag','outputFilenameBox'),'String',[fp fn]);
handles.outputFilename = [fp fn];
guidata(hObject,handles);

% --- Executes on button press in registerButton.
function registerButton_Callback(hObject, eventdata, handles)
% hObject    handle to registerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DW_Elastix_Transform(handles.movingFile,handles.outputFilename,handles.transformParameters);