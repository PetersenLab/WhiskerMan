function varargout = help_wm3(varargin)
% HELP_WM3 MATLAB code for help_wm3.fig
%      HELP_WM3, by itself, creates a new HELP_WM3 or raises the existing
%      singleton*.
%
%      H = HELP_WM3 returns the handle to a new HELP_WM3 or the handle to
%      the existing singleton*.
%
%      HELP_WM3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELP_WM3.M with the given input arguments.
%
%      HELP_WM3('Property','Value',...) creates a new HELP_WM3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before help_wm3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to help_wm3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help help_wm3

% Last Modified by GUIDE v2.5 15-Feb-2016 13:56:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @help_wm3_OpeningFcn, ...
                   'gui_OutputFcn',  @help_wm3_OutputFcn, ...
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


% --- Executes just before help_wm3 is made visible.
function help_wm3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to help_wm3 (see VARARGIN)

% Choose default command line output for help_wm3
handles.output = hObject;


% Extract input pmtrs 
if nargin == 4
    text = varargin{1};
else
    error(sprintf('Incorrect number of inputs (%d)', nargin))
end

set(handles.text_help,'String',text)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes help_wm3 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = help_wm3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1);

% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1)
