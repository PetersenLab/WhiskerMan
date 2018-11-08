function varargout = setroi_gui_wm4b(varargin)
% SETROI_GUI_WM4B MATLAB code for setroi_gui_wm4b.fig
%      SETROI_GUI_WM4B, by itself, creates a new SETROI_GUI_WM4B or raises the existing
%      singleton*.
%
%      H = SETROI_GUI_WM4B returns the handle to a new SETROI_GUI_WM4B or the handle to
%      the existing singleton*.
%
%      SETROI_GUI_WM4B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETROI_GUI_WM4B.M with the given input arguments.
%
%      SETROI_GUI_WM4B('Property','Value',...) creates a new SETROI_GUI_WM4B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before setroi_gui_wm4b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to setroi_gui_wm4b_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setroi_gui_wm4b

% Last Modified by GUIDE v2.5 06-Apr-2016 16:52:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @setroi_gui_wm4b_OpeningFcn, ...
                   'gui_OutputFcn',  @setroi_gui_wm4b_OutputFcn, ...
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


% --- Executes just before setroi_gui_wm4b is made visible.
function setroi_gui_wm4b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to setroi_gui_wm4b (see VARARGIN)

% Choose default command line output for setroi_gui_wm4b
handles.output = hObject;

% Extract input pmtrs (video frame and current ROIs)
if nargin == 5
    tmp = varargin{1};
    handles.frame.h = tmp.h.frame;
    handles.frame.v = tmp.v.frame;
    clear tmp
    handles.roi = varargin{2};
else
    error(sprintf('Incorrect number of inputs (%d)', nargin))
end

handles.linelen = 20;

axes(handles.hframe)
image(handles.frame.h)
hold on
x = handles.roi.h(1:2);
y = handles.roi.h(3:4);
plot([x(1) x(1) x(2) x(2) x(1)],[y(1) y(2) y(2) y(1) y(1)],'r--')

axes(handles.vframe)
image(handles.frame.v)
hold on
x = handles.roi.v(1:2);
y = handles.roi.v(3:4);
plot([x(1) x(1) x(2) x(2) x(1)],[y(1) y(2) y(2) y(1) y(1)],'b--')
clear x y

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes setroi_gui_wm4b wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = setroi_gui_wm4b_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

varargout{1} = handles.roi;
delete(handles.figure1);

% --- Executes on button press in pushbutton_setroih.
function pushbutton_setroih_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setroih (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.hframe)
cla
image(handles.frame.h)
hold on

topleft = round(ginput(1));
handles.roi.h(1) = topleft(1);
handles.roi.h(3) = topleft(2);
plot(handles.roi.h(1),handles.roi.h(3),'r+')

bottomright = round(ginput(1));
handles.roi.h(2) = bottomright(1);
handles.roi.h(4) = bottomright(2);
plot(handles.roi.h(2),handles.roi.h(4),'m+')

x = handles.roi.h(1:2);
y = handles.roi.h(3:4);
plot([x(1) x(1) x(2) x(2) x(1)],[y(1) y(2) y(2) y(1) y(1)],'r--')
clear x y

guidata(hObject, handles);

% --- Executes on button press in pushbutton_setroiv.
function pushbutton_setroiv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setroiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.vframe)
cla
image(handles.frame.v)
hold on

topleft = round(ginput(1));
handles.roi.v(1) = topleft(1);
handles.roi.v(3) = topleft(2);
plot(handles.roi.v(1),handles.roi.v(3),'b+')

bottomright = round(ginput(1));
handles.roi.v(2) = bottomright(1);
handles.roi.v(4) = bottomright(2);
plot(handles.roi.v(2),handles.roi.v(4),'g+')

x = handles.roi.v(1:2);
y = handles.roi.v(3:4);
plot([x(1) x(1) x(2) x(2) x(1)],[y(1) y(2) y(2) y(1) y(1)],'b--')
clear x y

guidata(hObject, handles);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1)
