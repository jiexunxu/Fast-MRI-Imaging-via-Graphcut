function varargout = DICOMViewer(varargin)
% DICOMViewer M-file for DICOMViewer.fig
%      DICOMViewer, by itself, creates a new DICOMViewer or raises the existing
%      singleton*.
%
%      H = DICOMViewer returns the handle to a new DICOMViewer or the handle to
%      the existing singleton*.
%
%      DICOMViewer('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICOMViewer.M with the given input arguments.
%
%      DICOMViewer('Property','Value',...) creates a new DICOMViewer or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DICOMViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DICOMViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% original from Michael Wunder submission id#4172
% modified by Olivier Salvado june04, Case Western Reserve University
% to rename the files using the dicom headers.

% Edit the above text to modify the response to help DICOMFILES

% Last Modified by GUIDE v2.5 16-Jun-2004 09:49:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DICOMViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @DICOMViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before DICOMFiles is made visible.
function DICOMViewer_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for DICOMFiles
handles.output = hObject;
handles.hidden = [];
guidata(hObject, handles);
colormap gray(256)

% use push-button callback
if length(varargin) & ischar (varargin{1})
   handles.dfolder = varargin{1};
   SetFolder(handles);
   ListBox_Callback(hObject, eventdata, handles);
else
   newFolder_Callback(hObject, eventdata, handles)
end

% UIWAIT makes DICOMFiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = DICOMViewer_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in newFolder.
function newFolder_Callback(hObject, eventdata, handles)

P = fileparts(mfilename('fullpath'));
% nfolder=uigetdir(P,'Select DICOM Directory');
nfolder=uigetdir('','Select DICOM Directory');
if ~ischar(nfolder)
    disp('no valid Directory selected.')
    return;
end
handles.dfolder=nfolder;
guidata(hObject, handles);
SetFolder(handles);
ListBox_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetFolder (handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dfiles=dir(handles.dfolder);
dfiles=dfiles(3:end);                   % avoid . and ..
nfiles=length(dfiles);
if nfiles<1
    disp('no files availabel.')
    return;
end
set(handles.ListBox,'String',char(dfiles.name),'value',1);
s = [num2str(nfiles) ' files in: ' handles.dfolder];
set(handles.NofFiles,'String', s);
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function ListBox_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in ListBox.
function ListBox_Callback(hObject, eventdata, handles)

fname = get(handles.ListBox,'String');
fname = fname(get(handles.ListBox,'value'),:);
try
   metadata = dicominfo([handles.dfolder '\' fname]);
catch
   disp ('apparently not a DICOM file');
   return
end
img      = dicomread([handles.dfolder '\' fname]);
imagesc(img);axis image
axis off

ch = get(handles.HeaderList, 'value');
fields=char(fieldnames(metadata));
len = setdiff (1:size(fields,1), handles.hidden);
id=0;
for k=len,
    estr=eval(['metadata.' fields(k,:)]);
    if ischar(estr)
        str=[fields(k,:) ' : ' estr];
    elseif isnumeric(estr)
        str=[fields(k,:) ' : ' num2str(estr(1:min(3,end))')];
    else
        str=[fields(k,:) ' : ...'];
    end
    id = id+1;
    cstr{id}=sprintf('%3d %s',k,str);
end
set(handles.HeaderList,'Value',ch);
set(handles.HeaderList,'String',cstr);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function HeaderList_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in cmdHide.
function cmdHide_Callback(hObject, eventdata, handles)
st = get (handles.HeaderList, 'string');
hide = get(handles.HeaderList, 'value');
if length(hide)==length(st)
   disp ('WARNING: at least one field must be shown');
   return
end
hidev=[];
for id=hide
   hidev = [hidev str2num(st{id}(1:3))];
end
handles.hidden = union (handles.hidden, hidev);
set(handles.HeaderList,'Value',1);
guidata (hObject, handles);
ListBox_Callback(hObject, eventdata, handles);

% --- Executes on button press in cmdShowAll.
function cmdShowAll_Callback(hObject, eventdata, handles)

handles.hidden = [];
guidata (hObject, handles);
ListBox_Callback(hObject, eventdata, handles);


% --- Executes on button press in ButtonRename.
function ButtonRename_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = get(handles.ListBox,'String');
h = waitbar(0,'Copying files....');
Nfiles = length(fname);
pid_raw_old = '';
for k=1:Nfiles,
    waitbar(k/Nfiles);
    try
        fnamek = fname(k,:);
        source = [handles.dfolder '\' fnamek];
        
        metadata = dicominfo(source);
        
        % --- get the fileds for the new file name
        pid_raw = char(metadata.PatientID);
        
        % check the pid
        if (length(pid_raw)>10) && (~strcmp(pid_raw,pid_raw_old)),
            prompt={'Enter the patient ID'};
            name='Wrong patient ID';
            numlines=1;
            defaultanswer={pid_raw};
            pid = char(inputdlg(prompt,name,numlines,defaultanswer));
        elseif (length(pid_raw)<=10)
            pid = pid_raw;
        end
        pid_raw_old = pid_raw;

        
%         seq = metadata.SequenceName;
        seq = metadata.SeriesDescription;
        ser = metadata.SeriesNumber;
        ima = num2str(metadata.InstanceNumber);
        TE = round(metadata.EchoTime);
        TR = round(metadata.RepetitionTime);
        b = '000';
        b(end-length(ima)+1:end) = ima;
        seq(isspace(seq)) = '_';

        % --- construct the new name
        newname = ['pid' pid '-seq' seq '-image' num2str(b) '-TE' num2str(TE)  '-TR' num2str(TR) '.dcm'];
        
        % --- create the a new directory with the series #
        wiam = cd;
        cd(handles.dfolder)
        dirname = ['PatientID' pid '-series' num2str(ser)];
        if exist([handles.dfolder '\' dirname],'dir')~=7,
            mkdir(dirname)
        end
        destination = [(handles.dfolder) '\' dirname '\' newname];
        
        % -- check if the file name exist, and add x's
        if exist(destination,'file'),
            while exist(destination,'file'),
                destination = [destination(1:end-4) 'x.dcm'];
            end
        end
            
        % --- copy the file
        [status,mess] = copyfile(source,destination,'f');
        
    catch
        disp ('apparently not a DICOM file, trying next one...');
        lasterr
    end
    
end

close(h)
