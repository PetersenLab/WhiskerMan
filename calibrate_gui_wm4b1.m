function varargout = calibrate_gui_wm4b1(varargin)
% CALIBRATE_GUI_WM4B1 MATLAB code for calibrate_gui_wm4b1.fig
%      CALIBRATE_GUI_WM4B1, by itself, creates a new CALIBRATE_GUI_WM4B1 or raises the existing
%      singleton*.
%
%      H = CALIBRATE_GUI_WM4B1 returns the handle to a new CALIBRATE_GUI_WM4B1 or the handle to
%      the existing singleton*.
%
%      CALIBRATE_GUI_WM4B1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATE_GUI_WM4B1.M with the given input arguments.
%
%      CALIBRATE_GUI_WM4B1('Property','Value',...) creates a new CALIBRATE_GUI_WM4B1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibrate_gui_wm4b1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibrate_gui_wm4b1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibrate_gui_wm4b1

% Last Modified by GUIDE v2.5 19-Aug-2016 10:29:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibrate_gui_wm4b1_OpeningFcn, ...
                   'gui_OutputFcn',  @calibrate_gui_wm4b1_OutputFcn, ...
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


% --- Executes just before calibrate_gui_wm4b1 is made visible.
function calibrate_gui_wm4b1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibrate_gui_wm4b1 (see VARARGIN)

% Choose default command line output for calibrate_gui_wm4b1
handles.output = hObject;


% Extract input pmtrs (video frame and current ROIs)
if nargin == 7
    handles.roi = varargin{1};
    tmp = varargin{2};
    handles.data_video.h.fname = tmp.h;
    handles.data_video.v.fname = tmp.v;
    tmp = varargin{3};
    handles.data_video.h.frame = tmp.h;
    handles.data_video.v.frame = tmp.v;
    tmp = varargin{4};
    handles.dir.h = tmp.h;
    handles.dir.v = tmp.v;
    clear tmp
else
    error(sprintf('Incorrect number of inputs (%d)', nargin))
end

ff = dir([handles.dir.h '*.dat']);
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose calibration video h';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_select_calib_video_h,'String',string_list);

ff = dir([handles.dir.v '*.dat']);
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose calibration video v';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_select_calib_video_v,'String',string_list);

axes(handles.axes_frame1)
cla
% xlim([1 calib_video.header.width])
% ylim([1 calib_video.header.height])
set(gca,'ydir','reverse')
imagesc(double(handles.data_video.h.frame(:,:,1)));
colormap gray
title(sprintf('%s',handles.data_video.h.fname))
axis image

axes(handles.axes_frame2)
cla
set(gca,'ydir','reverse')
imagesc(double(handles.data_video.v.frame(:,:,1)));
colormap gray
title(sprintf('%s',handles.data_video.v.fname))
axis image

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibrate_gui_wm4b1 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calibrate_gui_wm4b1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
if isfield(handles,'calib')
    varargout{1} = handles.calib;
else
    varargout{1} = handles.output;
end

delete(handles.figure1);

% --- Executes on selection change in popupmenu_select_calib_video_h.
function popupmenu_select_calib_video_h_Callback(hObject, eventdata, handles)

% Extract the selected file name:
val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    calib_video.h.fname = string_list{val};   % dat file.
else
    error('');
end
clear val string_list

fname = [handles.dir.h  calib_video.h.fname];
handles.calib_video.h = select_calib_video(fname,handles.axes_frame3);
handles.calib_video.h.fname = calib_video.h.fname;
% keyboard
clear calib_video fname

guidata(hObject, handles);

% --- Executes on selection change in popupmenu_select_calib_video_v.
function popupmenu_select_calib_video_v_Callback(hObject, eventdata, handles)

% Extract the selected file name:
val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    calib_video.v.fname = string_list{val};   % dat file.
else
    error('')
end
clear val string_list

fname = [handles.dir.v  calib_video.v.fname];
handles.calib_video.v = select_calib_video(fname,handles.axes_frame4);
handles.calib_video.v.fname = calib_video.v.fname;
clear calib_video fname

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_select_calib_video_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_calib_video_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ff = dir([handles.dir.v '\*.dat']);
% 
% nr_files=size(ff,1);
% string_list=cell(nr_files+1,1);
% string_list{1}='Choose calibration video v';
% for i=1:nr_files
%     string_list{i+1}=ff(i).name;
% end
% set(hObject,'String',string_list);

function calib_video = select_calib_video(fname,haxes)

% Open calib video file:
calib_video.fid = fopen(fname,'r');
calib_video.header = read_mikrotron_datfile_header(calib_video.fid);
calib_video.nframes = calib_video.header.nframes;
% video.width.raw = video.header.width;
% video.height.raw = video.header.height;
calib_video.offset = 8192;
% set file position to start of first frame
fseek(calib_video.fid,8192,-1);
% specify first frame
calib_video.startframe = calib_video.header.startframe+1;

% Load startframe:
offset = calib_video.header.imagesize * (calib_video.startframe-1) + calib_video.offset;
fseek(calib_video.fid,offset,-1);
clear offset
tmp = fread(calib_video.fid,calib_video.header.imagesize-24,'uint8=>uint8');
tmp = reshape([tmp; zeros(24,1)],calib_video.header.width,calib_video.header.height)';
calib_video.frame = uint8(zeros(calib_video.header.height,calib_video.header.width,3));
calib_video.frame(:,:,1) = tmp;
calib_video.frame(:,:,2) = tmp;
calib_video.frame(:,:,3) = tmp;
clear tmp

% Plot startframe:
axes(haxes)
cla
% xlim([1 handles.calib_video.header.width])
% ylim([1 handles.calib_video.header.height])
set(gca,'ydir','reverse')
imagesc(double(calib_video.frame(:,:,1)));
colormap gray
title(sprintf('%s',fname))
axis image

% 190816
% close the file, since it's only used here 
fclose(calib_video.fid);

% --- Executes during object creation, after setting all properties.
function popupmenu_select_calib_video_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_calib_video_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ff = dir([handles.dir.h '\*.dat']);
% 
% nr_files=size(ff,1);
% string_list=cell(nr_files+1,1);
% string_list{1}='Choose calibration video h';
% for i=1:nr_files
%     string_list{i+1}=ff(i).name;
% end
% set(hObject,'String',string_list);

% --- Executes on button press in pushbutton_check_alignment.
function pushbutton_check_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_check_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
data.h = handles.data_video.h.frame(:,:,1);
calib.h = handles.calib_video.h.frame(:,:,1);
data.v = handles.data_video.v.frame(:,:,1);
calib.v = handles.calib_video.v.frame(:,:,1);

if ~all(size(data.h)==size(calib.h)) || ~all(size(data.v)==size(calib.v))
    title('Image dimensions do not match')
    disp('Data image size:')
    size(data.h),size(data.v)
    disp('Calib image size:')
    size(calib.h),size(calib.v)
    beep
else
    subplot 121
    imagesc(data.h - calib.h), title('horizontal')
    subplot 122
    imagesc(data.v - calib.v), title('vertical')
end


% --- Executes on selection change in popupmenu_select_calib_file.
function popupmenu_select_calib_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_calib_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_calib_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_calib_file

val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    handles.calib_file = string_list{val};   % dat file.
else
    error('')
end
clear val string_list

% Apply ROIs, reflect the V view and plot them:
calib_frame.h = handles.calib_video_frame(handles.roi.h(3):handles.roi.h(4),handles.roi.h(1):handles.roi.h(2),:);
calib_frame.v = handles.calib_video_frame(handles.roi.v(3):handles.roi.v(4),handles.roi.v(1):handles.roi.v(2),:);
tmp = calib_frame.v(:,:,1)';
calib_frame.v = zeros(size(tmp,1),size(tmp,2),3);
calib_frame.v(:,:,1) = tmp;
calib_frame.v(:,:,2) = tmp;
calib_frame.v(:,:,3) = tmp;
handles.calib_frame = calib_frame;
clear tmp

axes(handles.axes_frame4)
cla
imagesc(double(calib_frame.h(:,:,1)));
axis image off
colormap gray
xlim([1 handles.calib_video.header.width])
ylim([1 handles.calib_video.header.height])
set(gca,'ydir','reverse')

axes(handles.axes_frame6)
cla
imagesc(double(calib_frame.v(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
ylim([1 handles.calib_video.header.width])
xlim([1 handles.calib_video.header.height])

data_frame.h = handles.data_video_frame(handles.roi.h(3):handles.roi.h(4),handles.roi.h(1):handles.roi.h(2),:);
data_frame.v = handles.data_video_frame(handles.roi.v(3):handles.roi.v(4),handles.roi.v(1):handles.roi.v(2) ,:);
tmp = data_frame.v(:,:,1)';
data_frame.v = zeros(size(tmp,1),size(tmp,2),3);
data_frame.v(:,:,1) = tmp;
data_frame.v(:,:,2) = tmp;
data_frame.v(:,:,3) = tmp;
clear tmp

axes(handles.axes_frame3)
cla
imagesc(double(data_frame.h(:,:,1)));
axis image off
colormap gray
xlim([1 handles.calib_video.header.width])
ylim([1 handles.calib_video.header.height])
set(gca,'ydir','reverse')

axes(handles.axes_frame5)
cla
imagesc(double(data_frame.v(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
ylim([1 handles.calib_video.header.width])
xlim([1 handles.calib_video.header.height])

% Now use 'C' to estimate the affine transform from real-world 3D coords to
% the H and V (orthographic) projection coordinates
load(handles.calib_file,'-mat','C')
if ~exist('C','var')
    error('Calib file must contain matrix ''C''.')
end
% C is in raw image coords.  Translate xy and XY components to reference frame of
% the xy ROI and vw components to reference frame of the vw ROI:
origin_XY = handles.roi.h([1,3]);
origin_xy = handles.roi.h([1,3]);   
origin_vw = handles.roi.v([1,3]);
C(:,1) = C(:,1) - origin_xy(1);
C(:,5) = C(:,5) - origin_XY(1);
C(:,2) = C(:,2) - origin_xy(2);
C(:,6) = C(:,6) - origin_XY(2);
C(:,3) = C(:,3) - origin_vw(1);
C(:,4) = C(:,4) - origin_vw(2);
% Rotate v-view calib coords:
tmp = C(:,[3,4]);
C(:,[3,4]) = tmp(:,[2,1]);
clear tmp

% Do the regression for 'v':
n = size(C,1);
v = C(:,3);
w = C(:,4);
P = C(:,5:7);
b = regress(v,[P ones(n,1)]);
calib.mv = b(1:3)';
calib.ov = b(4);
clear b
b = regress(w,[P ones(n,1)]);
calib.mw = b(1:3)';
calib.ow = b(4);
clear b
calib.matrix = [calib.mv;calib.mw];
calib.vector = [calib.ov;calib.ow];

% Save the results to 'init' file:
% roi = handles.roi;
% handles.init_file = [handles.calib_video_fname(1:end-3) 'init'];
% if exist(handles.init_file)
%     save(handles.init_file,'-mat','calib','roi','-append')
% else
%     save(handles.init_file,'-mat','calib','roi')
% end

%% Test regression consistency:
axes(handles.axes_frame6)
hold on
%assume frame 1 is plotted
frameidx = 1;
P1 = P(frameidx,:);
v1 = v(frameidx,:);
w1 = w(frameidx,:);
P2 = P(size(P,1)/2+frameidx,:);
v2 = v(size(P,1)/2+frameidx,:);
w2 = w(size(P,1)/2+frameidx,:);
plot(v1,w1,'m.',P1*calib.mv'+calib.ov,P1*calib.mw'+calib.ow,'ro')
plot(v2,w2,'m.',P2*calib.mv'+calib.ov,P2*calib.mw'+calib.ow,'ro')
clear P1 P2

handles.calib = calib;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_select_calib_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_calib_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.calib');

nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Select .calib file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);

% --- Executes on button press in pushbutton2.
function pushbutton_check_correspondence_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% user selects a point in the H view
figure
subplot 121
imagesc(double(handles.calib_frame.h(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
title('Horizontal view: Choose a point')
[x,y] = ginput(1);

% Now calculate the corresponding line in the v view
% For theory, see notebook entry 110216
% coefficents of linear polynomial in V view:
p.v(1) = handles.calib.mv(3);
p.v(2) = handles.calib.mv(1)*x + handles.calib.mv(2)*y + handles.calib.ov;
p.w(1) = handles.calib.mw(3);
p.w(2) = handles.calib.mw(1)*x + handles.calib.mw(2)*y + handles.calib.ow;
zmax = (1 - p.w(2))/p.w(1);
% zmin = (size(handles.calib_video_frame.v,1) - p.w(2))/p.w(1);
zmin = (handles.roi.v(4) - p.w(2))/p.w(1);
z = linspace(zmin,zmax,100);

% Plot the outcome
subplot 122
imagesc(double(handles.calib_frame.v(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
% ylim([1 handles.calib_video.header.width])
% xlim([1 handles.calib_video.header.height])
hold on
plot(polyval(p.v,z),polyval(p.w,z),'w')
title('Vertical view')


% --- Executes on button press in pushbutton_help_select_calib_video.
function pushbutton_help_select_calib_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_select_calib_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('First step is to choose a calibration video.  This *must* be obtained under imaging conditions identical to that of the data video selected in WhiskerMan.');

% --- Executes on button press in pushbutton_help_check_alignment.
function pushbutton_help_check_alignment_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_check_alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('Subtract frame of calibration video from frame of data video.  Allows you to check for subtle misalignments between the two.');

% --- Executes on button press in pushbutton_help_select_calib_file.
function pushbutton_help_select_calib_file_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_select_calib_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('Use calibration data to compute the orthographic transform from 3D real-world coordinates to the horizontal and vertical projections, by multiple regression.  Select a .calib file (a -Mat file).  This *must* correspond to the chosen calibration video.  This file must contain a matrix C with the following structure.  Each row specifies the location of a given point, by reference to 3 frames of reference: (1) the world (3 values); (2) horizontal projection (2 values) and (3) vertical projection (2 values).  These values are currently obtained from Mat Evan''s program texture_calibrate.');

% --- Executes on button press in pushbutton_help_check_correspondence.
function pushbutton_help_check_correspondence_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_check_correspondence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('To check calibration, assume that calibration frame has a given 3D point that is uniquely identifiable in both horizontal and vertical views.  This might, for example, be the end of a cut whisker, a water droplet etc.  When prompted, select a point in the horizontal view.  Using the calibration affine transform, the program computes the line in the vertical view corresponding to this point.  If calibration has worked correctly, this line will pass through the image of the target point in the vertical view.')

% --- Executes on button press in pushbutton_help_texture_calibrate.
function pushbutton_help_texture_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_texture_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('Use Mat Evan''s ''texture_calibrate'' function to extract calibration points from the currently selected calibration video.  The function saves the coordinates of these points as a matrix ''C'' with the format detailed in the help for ''Select calib file''.  The ''C'' matrix is saved in a file with extension ''.calib''.  The function uses the options ''pause'' and ''display'' which can be set using the corresponding checkboxes.')

% --- Executes on button press in pushbutton_help_texture_calibrate_options.
function pushbutton_help_texture_calibrate_options_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help_texture_calibrate_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

help_wm3('Options for ''texture_calibrate''.  Default is for the program to zip through all frames without either displaying. If ''display'' is checked, each frame is displayed.  If ''pause'' is checked, algorithm pauses after each frame')

% --- Executes on button press in pushbutton_texture_calibrate.
function pushbutton_texture_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_texture_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prepare input arguments for 'texture_calibrate':
% keyboard
fname.h = [handles.dir.h handles.calib_video.h.fname];
fname.v = [handles.dir.v handles.calib_video.v.fname];
framenums = round(linspace(1,handles.calib_video.h.nframes,100));
% keyboard
roi1 = handles.texture_calibrate_rois.h;
roi2 = handles.texture_calibrate_rois.v;
% roi1 = handles.roi.h;
% roi2 = handles.roi.v;
pause = get(handles.checkbox_texture_calibrate_pause,'value');
draw = get(handles.checkbox_texture_calibrate_display,'value');

figure
[C,image] = texture_calibrate_rsp_wm4b1(fname,framenums,draw,pause,roi1,roi2);
% C has 7 columns: (x,y,v,w,X,Y,Z)
% each row of C is the coordinates of a given point in H view (x,y), V view
% (v,w) and real-world 3D (X,Y,Z).
% x,y,v,w are pixel coordinates in the frame of the raw video image.
% Z comes from prior knowledge of the calibration object geometry, built in to the function.

% handles.calib_file = [fname(1:end-3) 'calib'];
% save(handles.calib_file,'-mat','C')


% % update file list in the calib select button:
% ff = dir('*.calib');
% nr_files=size(ff,1);
% string_list=cell(nr_files+1,1);
% string_list{1}='Select .calib file';
% for i=1:nr_files
%     string_list{i+1}=ff(i).name;
% end
% set(handles.popupmenu_select_calib_file,'String',string_list);

handles.texture_calibrate_results.C = C;
handles.texture_calibrate_results.image = image;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in checkbox_texture_calibrate_pause.
function checkbox_texture_calibrate_pause_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_texture_calibrate_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_texture_calibrate_pause


% --- Executes on button press in checkbox_texture_calibrate_display.
function checkbox_texture_calibrate_display_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_texture_calibrate_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_texture_calibrate_display

% --- Executes on button press in pushbutton_set_ROIs_for_texture_calibrate.
function pushbutton_set_ROIs_for_texture_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_ROIs_for_texture_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'texture_calibrate_rois')
    roi_default = handles.texture_calibrate_rois;
else
    roi_default = handles.roi;
end
handles.texture_calibrate_rois = setroi_gui_wm4b(handles.calib_video,roi_default);
% keyboard
clear roi_default
guidata(hObject, handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

calib = handles.calib;
C = handles.texture_calibrate_results;
calib_file = [handles.calib_video.h.fname(1:end-3) 'calib'];
save(calib_file,'calib','C')
clear calib C


uiresume(handles.figure1)


% --- Executes on button press in pushbutton_derive_calibration_transform.
function pushbutton_derive_calibration_transform_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_derive_calibration_transform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.texture_calibrate_results.C
if ~isfield(handles.texture_calibrate_results,'C')
    error('No ''C'' matrix.  Perhaps you have not ''Run texture calibrate''?')
else
    C = handles.texture_calibrate_results.C;
end

% Apply ROIs and plot them:
calib_frame.h = handles.calib_video.h.frame(handles.roi.h(3):handles.roi.h(4),handles.roi.h(1):handles.roi.h(2),:);
calib_frame.v = handles.calib_video.v.frame(handles.roi.v(3):handles.roi.v(4),handles.roi.v(1):handles.roi.v(2),:);
% tmp = calib_frame.v(:,:,1)';
% calib_frame.v = zeros(size(tmp,1),size(tmp,2),3);
% calib_frame.v(:,:,1) = tmp;
% calib_frame.v(:,:,2) = tmp;
% calib_frame.v(:,:,3) = tmp;
handles.calib_frame = calib_frame;
clear tmp

axes(handles.axes_frame4)
cla
imagesc(double(calib_frame.h(:,:,1)));
axis image off
colormap gray
xlim([1 handles.calib_video.h.header.width])
ylim([1 handles.calib_video.h.header.height])
set(gca,'ydir','reverse')

axes(handles.axes_frame6)
cla
imagesc(double(calib_frame.v(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
ylim([1 handles.calib_video.v.header.width])
xlim([1 handles.calib_video.v.header.height])

data_frame.h = handles.data_video.h.frame(handles.roi.h(3):handles.roi.h(4),handles.roi.h(1):handles.roi.h(2),:);
data_frame.v = handles.data_video.v.frame(handles.roi.v(3):handles.roi.v(4),handles.roi.v(1):handles.roi.v(2) ,:);
% tmp = data_frame.v(:,:,1)';
% data_frame.v = zeros(size(tmp,1),size(tmp,2),3);
% data_frame.v(:,:,1) = tmp;
% data_frame.v(:,:,2) = tmp;
% data_frame.v(:,:,3) = tmp;
% clear tmp

axes(handles.axes_frame3)
cla
imagesc(double(data_frame.h(:,:,1)));
axis image off
colormap gray
xlim([1 handles.calib_video.h.header.width])
ylim([1 handles.calib_video.h.header.height])
set(gca,'ydir','reverse')

axes(handles.axes_frame5)
cla
imagesc(double(data_frame.v(:,:,1)));
set(gca,'ydir','reverse')
axis image off
colormap gray
ylim([1 handles.calib_video.v.header.width])
xlim([1 handles.calib_video.v.header.height])

% Now use 'C' to estimate the affine transform from real-world 3D coords to
% the H and V (orthographic) projection coordinates
% load(handles.calib_file,'-mat','C')
% if ~exist('C','var')
%     error('Calib file must contain matrix ''C''.')
% end
% C is in raw image coords.  Translate xy and XY components to reference frame of
% the xy ROI and vw components to reference frame of the vw ROI:
origin_XY = handles.roi.h([1,3]);
origin_xy = handles.roi.h([1,3]); 
origin_vw = handles.roi.v([1,3]);

save('C.mat','C')
C(:,1) = C(:,1) - origin_xy(1);
C(:,5) = C(:,5) - origin_XY(1);
C(:,2) = C(:,2) - origin_xy(2);
C(:,6) = C(:,6) - origin_XY(2);
C(:,3) = C(:,3) - origin_vw(1);
C(:,4) = C(:,4) - origin_vw(2);
% % Rotate v-view calib coords:
% tmp = C(:,[3,4]);
% C(:,[3,4]) = tmp(:,[2,1]);
% clear tmp

% Do the regression for 'v':
n = size(C,1);
P = C(:,5:7);   % the 3D "real world" coords
v = C(:,3);     % projection of P onto V image plane
w = C(:,4);     % ditto
b = regress(v,[P ones(n,1)]);
calib.mv = b(1:3)';
calib.ov = b(4);
clear b
b = regress(w,[P ones(n,1)]);
calib.mw = b(1:3)';
calib.ow = b(4);
clear b
calib.matrix = [calib.mv;calib.mw];
calib.vector = [calib.ov;calib.ow];

%% Test regression consistency:
axes(handles.axes_frame6)
hold on
%assume frame 1 is plotted
frameidx = 1;
P1 = P(frameidx,:);
v1 = v(frameidx,:);
w1 = w(frameidx,:);
P2 = P(size(P,1)/2+frameidx,:);
v2 = v(size(P,1)/2+frameidx,:);
w2 = w(size(P,1)/2+frameidx,:);
plot(v1,w1,'m.',P1*calib.mv'+calib.ov,P1*calib.mw'+calib.ow,'ro')
plot(v2,w2,'m.',P2*calib.mv'+calib.ov,P2*calib.mw'+calib.ow,'ro')
clear P1 P2

%% Finish up:

handles.calib = calib;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_quit__without_save.
function pushbutton_quit__without_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_quit__without_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1)
