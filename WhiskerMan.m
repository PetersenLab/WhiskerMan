function varargout = WhiskerMan(varargin)
% WHIKERMAN2L M-file for WhikerMan2l.fig
%      WHIKERMAN2L, by itself, creates a new WHIKERMAN2L or raises the
%    .  existing
%      singleton*.
%
%      H = WHIKERMAN2L returns the handle to a new WHIKERMAN2L or the
%      handle to
%      the existing singleton*.
%4
%      WHIKERMAN2L('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WHIKERMAN2L.M with the given input arguments.
%
%      WHIKERMAN2L('Property','Value',...) creates a new WHIKERMAN2L or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before WhikerMan2l_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to WhikerMan2l_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WhikerMan2l

% Last Modified by GUIDE v2.5 08-Nov-2018 11:31:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WhikerMan2l_OpeningFcn, ...
    'gui_OutputFcn',  @WhikerMan2l_OutputFcn, ...
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


% --- Executes just before WhikerMan2l is made visible.
function WhikerMan2l_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WhikerMan2l (see VARARGIN)

% Choose default command line output for WhikerMan2l
handles.output = hObject;



% default pmtrs
handles.red = [1 0 0]; handles.green = [0 1 0];

handles.trpmtrs.sigma_prior = .1;    % ditto
set(handles.edit_sigma_prior,'String',handles.trpmtrs.sigma_prior);

handles.trpmtrs.sigma2_prior = .001;    % ditto
set(handles.edit_sigma2_prior,'String',handles.trpmtrs.sigma2_prior);


handles.default_energy_threshold_subtract_image_mean = -20;
handles.default_energy_threshold_do_not_subtract_image_mean = 80;

handles.frame_interval = 1;
set(handles.edit_frame_interval,'String',handles.frame_interval)

set(handles.edit_snout_sigma,'String',20)

handles.continuous_tracking = false;
set(handles.radiobutton_continuous_tracking,'value',handles.continuous_tracking)

handles.trpmtrs.subtract_image_mean = false;
%set(handles.radiobutton_subtract_image_mean,'value',handles.trpmtrs.subtract_image_mean)

% new in 4b2 (050916):
handles.trpmtrs.tracking_direction = 'fwd';
set(handles.pushbutton_tracking_direction,'String',handles.trpmtrs.tracking_direction)

%%%define options for the solver (warning in matlab 2017) Andrea
handles.options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'TolX',.05,'Display','off');

% if get(handles.radiobutton_subtract_image_mean,'Value')
%     handles.energy_threshold = handles.default_energy_threshold_subtract_image_mean;
% else
%     handles.energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
% end
% set(handles.edit_energy_threshold,'String',handles.energy_threshold);

% handles.trpmtrs.miumax = 0.2;
% set(handles.edit_miumax,'String',handles.trpmtrs.miumax);

% % default [x1,x2,y1,y2] coords where x1<x2, y1<y2:
% handles.roi_default.h = [300 700 1 320];
% handles.roi_default.v = [1 250 100 320];

% set(handles.checkbox_rotatevview,'value',1)

% handles.Nwhiskers = 1;
% set(handles.edit_Nwhiskers,'String',handles.Nwhiskers);

% User is not allowed to start tracking, until calibration is done:
set(handles.pushbutton_fit_snakes,'enable','off')

% set(handles.edit_playback_pause,'String',50)

% colour code for plotting different whiskers:
handles.colours = {'b.-','g.-','r.-','c.-','y.-','m.-','w.-','k.-'};
handles.lines = {'b-','g-','r-','c-','y-','m-','w-','k-'};
handles.lines_dotted = {'b:','g:','r:','c:','y:','m:','w:','k:'};
handles.points = {'b.','g.','r.','c.','y.','m.','w.','k.'};

handles.mastervideo_selected = false;   % mastervideo is the horizontal view

% new in v4b4:
set(handles.edit_snout_outliers,'String',3)

% new in v4b5:
%set(handles.edit_kinematics_test_position,'String','.5');

%new in v4c1:
handles.tracking2D=false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify folders in next 2 lines:
%set user
%user=andrea;
% switch user
%     case 'andrea'
        handles.dir.h =[pwd '\'];
        handles.dir.v =[pwd '\vertical_view\'];
%     case 'rasmus'
%         handles.dir.h ='C:\Users\mjcssrp\Dropbox\Whikerman4b5\ExampleH\';
%         handles.dir.v ='C:\Users\mjcssrp\Dropbox\Whikerman4b5\ExampleV\';
%     otherwise
%         error(['Unknown user ' user])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = handles.dir.h;
ff = dir([d '*.dat']);
clear d
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose H video';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_choose_hvideo,'String',string_list);
clear nr_files string_list i

d = handles.dir.v;
ff = dir([d '*.dat']);
clear d
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose V video';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_choose_vvideo,'String',string_list);
clear nr_files string_list i

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WhikerMan2l wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WhikerMan2l_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu_choose_hvideo.
function popupmenu_choose_hvideo_Callback(hObject, eventdata, handles)

% Until calib information has been provided, do not allow tracking:
if handles.tracking2D
set(handles.pushbutton_fit_snakes,'enable','on')
end
val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    handles.fname.h = string_list{val};   % avi file or dat file.
else
    error('')
end
clear val string_list

% get a FID and header information ('video' structure):
handles.video.h = initialise_new_video(handles.fname.h,handles,handles.dir.h);
% Set default ROIs:
handles.roi.h = [1,handles.video.h.width.raw,1,handles.video.h.height.raw];
handles.video.h.width.roi = handles.roi.h(2)-handles.roi.h(1)+1;
handles.video.h.height.roi = handles.roi.h(4)-handles.roi.h(3)+1;
handles.currentframe.h = handles.video.h.startframe;

handles.mastervideo_selected = true;

if handles.tracking2D
        % make the master frameidx table
%     if handles.video.h.nframes==handles.video.v.nframes
        nframes = handles.video.h.nframes;
        % straightforward
        frameidxs.h = modadd(1:nframes,handles.video.h.startframe-1,nframes);
%         frameidxs.v = modadd(1:nframes,handles.video.v.startframe-1,nframes);
        idx = find(frameidxs.h==1);
%         handles.framenums_vfromh = modadd(1:nframes,frameidxs.v(idx)-1,nframes);
        clear frameidxs idx
%     else
%         error('Extend the code!')
%     end
    
    % initialise parameters that apply at the level of the 2 views together
    handles.dt = 0.01;
    hsize = [3 3]; % default
    sigma = .05;
    handles.gaussian = fspecial('gaussian', hsize, sigma);
    clear hsize sigma
    % meanframe
    meanframe = init_subtract_image_mean(handles.trpmtrs.subtract_image_mean,handles.roi,handles.video);
    % 060416: is following line needed?
    fn = [handles.fname.h(1:end-3) 'meanframe'];    % NB H is the master
    save(fn,'meanframe')
    handles.meanframe = meanframe;
    clear fn meanframe
    % load and plot the first frame:
    titles.h = sprintf('Frame %d',handles.currentframe.h);
%     titles.v = sprintf('Frame %d',handles.currentframe.v);
    haxes.h = handles.viewh; 
%     haxes.v = handles.viewv; 
    handles.frame = load_and_plot_frames(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,haxes);
    clear haxes
    % set default tracking file
    %xx check for '.tr4' and modify according to next line...
    handles.trfname = [handles.fname.h(1:end-4) '.tr4_2D'];
    % If a 'tr' file for current video exists, load the data; else, initialise:
    handles = initialise_new_track(handles);
    handles.current_view = 1;
    handles.current_pt = 0;
    
    %xx add following line
    calib.matrix=zeros(2,3);
    calib.vector=zeros(2,1);
    calib.mv=zeros(1,3);
    calib.mw=zeros(1,3);
    calib.ov=0;
    calib.ow=0;
    handles.calib = calib;
    
    if isfield(handles,'folliclemask')    
    handles= rmfield(handles,'folliclemask');
    end
    if isfield(handles,'calib')
        plot_contours(handles.currentframe,handles)
    end
end


guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_choose_vvideo.
function popupmenu_choose_vvideo_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    handles.fname.v = string_list{val};   % avi file or dat file.
else
    error('')
end
clear val string_list

% get a FID and header information ('video' structure):
handles.video.v = initialise_new_video(handles.fname.v,handles,handles.dir.v);
handles.roi.v = [1,handles.video.v.width.raw,1,handles.video.v.height.raw];
handles.video.v.width.roi = handles.roi.v(2)-handles.roi.v(1)+1;
handles.video.v.height.roi = handles.roi.v(4)-handles.roi.v(3)+1;

handles.currentframe.v = handles.video.v.startframe;

% user should already have selected the "master" horizontal view video
if handles.mastervideo_selected
    % make the master frameidx table
    if handles.video.h.nframes==handles.video.v.nframes
        nframes = handles.video.h.nframes;
        % straightforward
        frameidxs.h = modadd(1:nframes,handles.video.h.startframe-1,nframes);
        frameidxs.v = modadd(1:nframes,handles.video.v.startframe-1,nframes);
        idx = find(frameidxs.h==1);
        handles.framenums_vfromh = modadd(1:nframes,frameidxs.v(idx)-1,nframes);
        clear frameidxs idx
    else
        error('Extend the code!')
    end
    
    % initialise parameters that apply at the level of the 2 views together
    handles.dt = 0.01;
    hsize = [3 3]; % default
    sigma = .05;
    handles.gaussian = fspecial('gaussian', hsize, sigma);
    clear hsize sigma
    % meanframe
    meanframe = init_subtract_image_mean(handles.trpmtrs.subtract_image_mean,handles.roi,handles.video);
    % 060416: is following line needed?
    fn = [handles.fname.h(1:end-3) 'meanframe'];    % NB H is the master
    save(fn,'meanframe')
    handles.meanframe = meanframe;
    clear fn meanframe
    % load and plot the first frame:
    titles.h = sprintf('Frame %d',handles.currentframe.h);
    titles.v = sprintf('Frame %d',handles.currentframe.v);
    haxes.h = handles.viewh; 
    haxes.v = handles.viewv; 
    handles.frame = load_and_plot_frames(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,haxes);
    clear haxes
    % set default tracking file
    handles.trfname = [handles.fname.h(1:end-4) '.tr4'];
    % If a 'tr' file for current video exists, load the data; else, initialise:
    handles = initialise_new_track(handles);
    handles.current_view = 1;
    handles.current_pt = 0;
    if isfield(handles,'folliclemask')    
    handles= rmfield(handles,'folliclemask');
    end
    if isfield(handles,'calib')
        plot_contours(handles.currentframe,handles)
    end
else
    error('Select the horizontal view (master) video first')
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_contours(framenum,handles)
% plot tracked contour for selected frame

frameidx = framenum.h;
dt = 0.05;
t = 0:dt:1;
for w = 1:handles.Nwhiskers
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    
    if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx)
        r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
        b3 = bezierval(r3,t);
        r2 = projection2(r3,handles.calib);
        b2 = projection2(b3,handles.calib);
        fp3 = handles.whisker(w).fp3_all(frameidx,:)';
        fp2 = projection2(fp3,handles.calib);
        axes(handles.viewh), hold on %#ok<LAXES>
        plot(b2(1,:,1),b2(2,:,1),handles.lines{w},r2(1,:,1),r2(2,:,1),handles.points{w},...
            fp2(1,1),fp2(2,1),'y.')

        if ~handles.tracking2D
            axes(handles.viewv), hold on %#ok<LAXES>
            plot(b2(1,:,2),b2(2,:,2),handles.lines{w},r2(1,:,2),r2(2,:,2),handles.points{w},...
            fp2(1,2),fp2(2,2),'y.')
        end
    end
end
clear w wmod r3 r2 b3 b2 fp3 fp2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b2 = projection(b3,calib)
% b3 ~ 3xN is a list of 3D coordinates
% b2 ~ projection of b3 onto the plane defined by 'calib'
N = size(b3,2);
b2 = calib.matrix*b3 + calib.vector*ones(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b2 = projection2(b3,calib)
% b3 ~ 3xNx2 is a list of 3D coordinates
% b2(:,:,2) ~ projection of b3 onto xy
% b2(:,:,2) ~ projection of b3 onto the plane defined by 'calib' (vw)
N = size(b3,2);
b2 = zeros(2,N,2);
b2(:,:,1) = b3(1:2,:);
if ~isempty(calib)
b2(:,:,2) = calib.matrix*b3 + calib.vector*ones(1,N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function handles = initialise_new_video(handles)
function video = initialise_new_video(fname,handles,dir)

% Fixed pmtrs (and therefore don't need to be loaded from tr file)

switch fname(end-2:end)
%    % To syncrhonise the cameras, need "startframe", so need more header info
%     case 'avi'
%         video.type = 'avi';
%         video.vObj = mmreader([dir fname]);
%         video.nframes = video.vObj.NumberOfFrames;
%         video.width.raw = video.vObj.Width;
%         video.height.raw = video.vObj.Height;
%         % specify first frame and last frame:
%         video.startframe = 1;
%         video.stopframe = video.nframes;
    case 'dat'
        video.type = 'dat';
        video.fid = fopen([dir fname],'r');
        video.header = read_mikrotron_datfile_header(video.fid);
        video.nframes = video.header.nframes;
        video.width.raw = video.header.width;
        video.height.raw = video.header.height;
        video.offset = 8192;
        % set file position to start of first frame
        fseek(video.fid,8192,-1);
        % specify first frame
        video.startframe = video.header.startframe+1;
        % new 070416 - replacement for previous line:
%         video.startframe = video.header.startframe+video.header.triggerframe;
        if video.startframe > 1
            video.stopframe = video.startframe - 1;
        else
            video.stopframe = video.nframes;
        end
        %video.header.startframe,video.header.triggerframe
    otherwise
        error(fprintf('Unhandled video file type: %s\n',fname(end-2:end)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function popupmenu_choose_hvideo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_hvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function popupmenu_choose_vvideo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_vvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% d = handles.dir.v;
% ff = [dir([d '\*.avi']);dir([d ' \*.dat'])];
% 
% nr_files=size(ff,1);
% string_list=cell(nr_files+1,1);
% string_list{1}='Choose V video';
% for i=1:nr_files
%     string_list{i+1}=ff(i).name;
% end
% set(hObject,'String',string_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_fit_snakes.
function pushbutton_fit_snakes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit_snakes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% is 'fnoutput needed'?
% [handles, fnoutput] = fit_snakes(handles);
[handles] = fit_snakes(handles);
clear fnoutput
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [handles,fnoutput] = fit_snakes(handles)
try
fnoutput=1;
switch handles.trpmtrs.tracking_direction   % v42b 050916
    case 'fwd'
        track_fwd = true;
    case 'bkwd'
        track_fwd = false;
    otherwise
        error(['Unrecognised handles.trpmtrs.tracking_direction ' handles.trpmtrs.tracking_direction])
end

% h video is master, v video slave.  So, eg, 'lastframe' is lastframe of h video. 
if track_fwd   % v42b 050916
    lastframe = modsubtract(handles.currentframe.h,1,handles.video.h.nframes);
else
    lastframe = modadd(handles.currentframe.h,1,handles.video.h.nframes);
end     
track = true;
firstframe = false(1,handles.Nwhiskers);

for w = 1:handles.Nwhiskers
    
    if handles.whisker(w).selected
        if (handles.currentframe.h==handles.video.h.startframe) || ~handles.whisker(w).tracked(lastframe)
            firstframe(w) = true;
        end
    end
end
clear w lastframe

doplot_none = 0;
doplot_light = 1;
doplot_full = 2;

while (track)
  
    track = get(handles.radiobutton_continuous_tracking,'Value');
    
    if track_fwd   % v42b 050916
        lastframe = modsubtract(handles.currentframe.h,1,handles.video.h.nframes);
        lastbutoneframe = modsubtract(handles.currentframe.h,2,handles.video.h.nframes);
    else
        lastframe = modadd(handles.currentframe.h,1,handles.video.h.nframes);
        lastbutoneframe = modadd(handles.currentframe.h,2,handles.video.h.nframes);
    end
    
    if any(firstframe) || (rem(handles.currentframe.h,handles.frame_interval)==0)
        doplot = doplot_light;
        % following 2 lines needed?
        %         axes(handles.viewv)
        %         cla reset
    else
        doplot = doplot_none;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocess horizontal:
    im.h.raw = double(handles.frame.h(:,:,1)) - handles.meanframe.h(:,:,1);
    % Smooth image a bit to facilitate whisker tracking
    im.h.s = imfilter(im.h.raw,handles.gaussian);
    [im.h.dx,im.h.dy] = gradient(im.h.s);
    
    if ~handles.tracking2D
    % Preprocess vertical:
    im.v.raw = double(handles.frame.v(:,:,1)) - handles.meanframe.v(:,:,1);
    im.v.s = imfilter(im.v.raw,handles.gaussian);
    [im.v.dx,im.v.dy] = gradient(im.v.s);
    else
        im.v=[];
    end
%     if doplot
%         fname = handles.fname(1:end-4);
%         idx = regexp(fname,'\_');
%         fname(idx) = '-';
%         clear idx
%         %         imagesc(im.s), colormap gray
%         %         title(sprintf('%s: Frame %d',fname,handles.currentframe))
%         %         hold on
% %         plot_image(im,handles)
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Whisker tracking
    
    D = handles.n_control_pts;
    
    % Snout outline detection:
    % if it's the firstframe, process within large Region Of Interest
    % (ROI).  else, pick an ROI centred around follicle position of
    % previous frame:
    
    roi.h.x = 50:size(handles.frame.h,2);%50
    roi.h.y = 50:size(handles.frame.h,1);%50

    if ~handles.tracking2D
    roi.v.x = 50:size(handles.frame.v,2);%50
    roi.v.y = 50:size(handles.frame.v,1);%50
    end

    
    [folliclemask{1}] = snout_segment(double(handles.frame.h(:,:,1)),roi.h,[],[],[],'horizontal',handles); 
    handles.folliclemask{handles.currentframe.h,1} = folliclemask{1};
    if ~handles.tracking2D
    [folliclemask{2}] = snout_segment(double(handles.frame.v(:,:,1)),roi.v,[],[],[],'vertical',handles);
    handles.folliclemask{handles.currentframe.h,2} = folliclemask{2};
    end
    clear folliclemask_old fpidx_old roi %theta
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(handles.radiobutton_plot_snout_contour,'value')
        axes(handles.viewh)
        plot(folliclemask{1}(1,:),folliclemask{1}(2,:),'y-')
        if ~handles.tracking2D
            axes(handles.viewv)
            plot(folliclemask{2}(1,:),folliclemask{2}(2,:),'y-')
        end
    end  
%     clear folliclemask

    goodautosolution = false(1,handles.Nwhiskers);
    if any(firstframe) && get(handles.radiobutton_auto_initialise,'value')
        % working on this in wm4a1...
        % set control pts automatically, using the tr data of the file
        % currently selected in the auto init popupmenu:
        menu = cellstr(get(handles.popupmenu_choose_auto_initialise_file,'String'));
        trfile = menu{get(handles.popupmenu_choose_auto_initialise_file,'Value')};
        load(trfile,'-mat','calib','whisker')
        handles.calib=calib;
        % Initialise whisker parameters from the tr file:
        if exist('whisker','var')
            handles.Nwhiskers = length(whisker);
            for w = 1:handles.Nwhiskers
                handles.whisker(w).label = whisker(w).label;
                handles.whisker(w).energy_threshold = whisker(w).energy_threshold;
                handles.whisker(w).selected = whisker(w).selected;
                handles.whisker(w).tracked=zeros(1,handles.video.h.nframes);
                
            end
            set(handles.text_Nwhiskers,'String',handles.Nwhiskers)
            handles.current_whisker = 1;
            update_current_whisker_display(handles.current_whisker, handles);
        else
            % parameters will be those set up by prior call to initialise_new_track
        end
        clear whisker 
        % initialise tracking parameters from the  tr file:
        if exist('trpmtrs','var')
            if ~isfield(trpmtrs,'tracking_direction')   % new in v4b2 050916
                % for backwards compatibility with files that might lack this pmtr
                trpmtrs.tracking_direction = handles.trpmtrs.tracking_direction;
            end
            handles.trpmtrs = trpmtrs;
            clear trpmtrs
            set(handles.edit_sigma_prior,'String',handles.trpmtrs.sigma_prior);
            if isfield(handles,'trpmtrs_sigma2_prior')
                set(handles.edit_sigma2_prior,'String',handles.trpmtrs_sigma2_prior);
            end
            if isfield(handles,'trpmtrs_snout_sigma')
                set(handles.edit_snout_sigma,'String',handles.trpmtrs_snout_sigma);
            end
            set(handles.radiobutton_subtract_image_mean,'value',handles.trpmtrs.subtract_image_mean)
            %set(handles.pushbutton_tracking_direction,'String',handles.trpmtrs_tracking_direction);
        end
        %update size of firstframe=Nwhiskers
        for w = 1:handles.Nwhiskers
            handles.whisker(w).selected=true;
                if (handles.currentframe.h==handles.video.h.startframe) || ~handles.whisker(w).tracked(lastframe)
                    firstframe(w) = true;
                end
        end
        [r3auto, goodautosolution] = auto_initialise(trfile,handles);
        %%if whikerman cannot find a good solution, then stops tracking and
        %%continues with a new video
        if any(goodautosolution)
            for w = 1:handles.Nwhiskers
                if ~goodautosolution(w)
                  handles.whisker(w).selected=false;
                  firstframe(w) = false;
                end
            end
        %%%%deselect whiskers with no good solution
        else  %%no good solution    
            fnoutput=2;
            disp('No initialisation')
            track=0;
            break
        end
        clear menu trfile
    end
   
    if any(firstframe)
        menu = cellstr(get(handles.popupmenu_initialise_from_file,'String'));
        trfile = menu{get(handles.popupmenu_initialise_from_file,'Value')};
        if strcmp(trfile(end-2:end),'tr4') || strcmp(trfile(end-5:end),'tr4_2D')
            plot_tracks_from_file(trfile,handles)
        end
        clear menu trfile
    end

    
    % loop over whiskers
    for w = 1:handles.Nwhiskers
        
        wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
        if handles.whisker(w).selected==false
            
            continue
        else
            % carry on
        end
        
%         titles.h = ''; titles.v = '';
%         plot_frame(handles.frame,handles.video,handles.meanframe,titles,handles)

        update_current_whisker_display(w, handles)
        r3all = handles.whisker(w).r3all;
        
        nattempts = 2;    % predict_contour needs work (2nd priority)       
        clear attempts
        for attempt = 1:nattempts
    
            switch attempt
                case 1
                    usepredictor = 1;
                case 2
                    usepredictor = 0;
                case 3
                    usepredictor = 2;
                    
                otherwise
                    error('Should never happen')
            end
            
            if firstframe(w)
                % no tracking solution for previous frame:
                % must set control pts. can be done manually or automatically,
                % depending on setting of the auto init radio button.
                r2old = zeros(2,D,2);    % x/y, CP, h/v
                r3old = zeros(3,D);      % x/y/z, CP
                if goodautosolution(w)
                    
                    r3old = r3auto{w};
                else
                    
                    % if either the autoinit failed or user has chosen not to use
                    % it, manually set the control pts.
                    % it's necessary that the first point is the one nearest the follicle and that the points thereafter move outwards
                    % monotonically along the whisker.
                    %
                    % first, set control points in the H view:
                    
                    axes(handles.viewh)
                    %                     title(sprintf('Initialise whisker %d (%s)',w,handles.whisker(w).label),'color',[1 0 0])
                    set(handles.currentframe_display,'String',sprintf('Initialise whisker %d (%s)',w,handles.whisker(w).label))
                    %                     defcol = get(handles.text_viewh,'backgroundcolor');
                    %                     set(handles.text_viewh,'backgroundcolor',[1 0 0])
                    for i = 1:D
                        [r2old(1,i,1),r2old(2,i,1)] = ginput(1);
                        plot(r2old(1,i,1),r2old(2,i,1),'wo')
                        r3old(1,i) = r2old(1,i,1);
                        r3old(2,i) = r2old(2,i,1);
                    end
                    %                     clear r2old
                    %                     set(handles.text_viewh,'backgroundcolor',defcol)
                    %                     title('')
                    %
                    %  second, use the H view CPs to constrain the V view ones:
                    axes(handles.viewv)
                    %                     title(sprintf('Initialise whisker %d (%s)',w,handles.whisker(w).label),'color',[1 0 0])
                    set(handles.currentframe_display,'String',sprintf('Initialise whisker %d (%s)',w,handles.whisker(w).label))
                    %                     defcol = get(handles.text_viewv,'backgroundcolor');
                    %                     set(handles.text_viewv,'backgroundcolor',[1 0 0])
                    clear ph
                    if ~handles.tracking2D
                    for i = 1:D
                        [pv,pw,z] = hview_point_2_vview_line(r2old(1,i,1),r2old(2,i,1),handles.calib,handles.roi);
                        % pv and pw are linear polyval coefficient vectors
                        % Plot the line that corresponds to the previously
                        % selected H view points:
                        ph(i) = plot(polyval(pv,z),polyval(pw,z),'w');
                        clear pv pw z
                        % Now user can use this to contrain initial points
                        % in V view:
                        [r2old(1,i,2),r2old(2,i,2)] = ginput(1);
                        plot(r2old(1,i,2),r2old(2,i,2),'wo')
                        % Since the H and V view points are corresponding,
                        % can solve the projection equation for 'z':
                        r3old(3,i) = z_from_xyvw(r2old(:,i,1),r2old(:,i,2),handles.calib.matrix, handles.calib.vector);
                    end
                    
                    delete(ph)
                    end
                    %                     clear r2old
                    %                     set(handles.text_viewv,'backgroundcolor',defcol)
                    %                     title('')
                    clear i defcol
                    usepredictor = 0; % first frame, so there is no previous frame to predict from
                    %                     titles.h = ''; titles.v = '';
                    %                     plot_frame(handles.frame,handles.video,handles.meanframe,titles,handles)
                    %                     axes(handles.viewh)
                    %                     plot(folliclemask{1}(1,:),folliclemask{1}(2,:),'y-')
                    %                     axes(handles.viewv)
                    %                     plot(folliclemask{2}(1,:),folliclemask{2}(2,:),'y-')
                    %                     for w2 = 1:w
                    %                         wmod = rem(w2,handles.Nwhiskers) + handles.Nwhiskers*(w2==handles.Nwhiskers);
                    %                         b3 = bezierval(squeeze(handles.whisker(w2).r3all(handles.currentframe,:,:)),0:.02:1);
                    %                         b2 = projection2(b3,handles.calib);
                    %                         axes(handles.viewh)
                    % %                         plot(r2old(1,:,1),r2old(2,:,1),'wo')
                    %                         plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod})
                    %                         axes(handles.viewv)
                    % %                         plot(r2old(1,:,2),r2old(2,:,2),'wo')
                    %                         plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod}),
                    %                     end
                    %                     clear w2 b3 b2
                end
                
            else
                r3old = reshape(r3all(lastframe,:,:),3,D);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Aim of next section of code is to set initial conditions for
            % the optimisation.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ddr3 = zeros(3,D);
            z = modsubtract(handles.currentframe.h,handles.video.h.startframe,handles.video.h.nframes);
            if z==handles.video.h.nframes
                z = 0;
            end
            if z >= 2
                if handles.whisker(w).tracked(lastframe) && handles.whisker(w).tracked(lastbutoneframe)
                    ddr3 = reshape(r3all(lastframe,:,:),3,D) - reshape(r3all(lastbutoneframe,:,:),3,D);
                end
            end
            
            % Various methods for setting the initial conditions:
            switch usepredictor
                case 1
                    % assume that contour is moving at constant velocity
                    % ('ddr3')
                    r3 = predict_contour_wm4(r3old,ddr3);
                case 0
                    % assume that contour is stationary
                    r3 = r3old;
                    %                 case 2
                    %                     % predict that contour is the closest one.  Search for nearest contour to rold
                    %                     r3 = find_contour_wm4(r3old,im.h, im.v, calib_matrix, calib_offset);
                    %                     ! the logic changes substantially in the move to 3D.
                    %                     ! Function needs to be generalised (2nd priority)
            end
            clear ddr3 %miumax
            
            if doplot==doplot_full
                r2 = projection2(r3,handles.calib);
                r2old = projection2(r3old,handles.calib);
                axes(handles.viewh)
                % plot solution from previous frame, if there is one:
                plot(r2old(1,:,1),r2old(2,:,1),'b+'),
                % plot initial conditions:
                plot(r2(1,:,1),r2(2,:,1),'c+')
                if ~handles.tracking2D
                axes(handles.viewv)
                plot(r2old(1,:,2),r2old(2,:,2),'b+'),
                plot(r2(1,:,2),r2(2,:,2),'c+')
                end
                clear r2 r2old
            end
            
            clear rold
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Optimise Bezier curve to best fit image and prior constraints
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Constrain P0 and P2 to move normal to contour at t=0 and t=1 respectively. 
            np0 = normal_plane(r3,0);   % np0 ~ column vectors that span the plane
            np1 = normal_plane(r3,1);
            
            % Optimisation parameters for 3D:
            % z(1:2) - components in plane normal to curve at P0=r3(:,1)
            % z(3:5) - P1=r3(:,2)
            % z(6:7) - components in plane normal to ourve at P2=r3(:,3)
            
            % Set initial values
            z0 = [0 0 r3(:,2)' 0 0];
            
            % Find 'z' that minimises the cost Bimfun(z)
            [z,Emin,exitflag,output] = fminunc(@(z) Bimfun_wm4(z, r3, np0, np1, im.h, im.v, handles.dt, ...
                handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib),z0, ....
                handles.options);
            clear z0 exitflag output
            
            % The solution:
            r3new = zeros(size(r3));
            r3new(:,1) = r3(:,1) + z(1)*np0(:,1) + z(2)*np0(:,2);
            r3new(:,2) = z(3:5)';
            r3new(:,3) = r3(:,3) + z(6)*np1(:,1) + z(7)*np1(:,2);
            clear z
            
            b3 = bezierval(r3new,0:.01:1);
            
            if doplot==doplot_full
                b2 = projection2(b3,handles.calib);
                r2 = projection2(r3,handles.calib);
                axes(handles.viewh)
                plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod},r2(1,:,1),r2(2,:,1), 'LineWidth',2,handles.points{wmod})
                if ~handles.tracking2D
                axes(handles.viewv),
                plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod},r2(1,:,2),r2(2,:,2), 'LineWidth',2,handles.points{wmod}),
                end
                clear b2 r2
            end
            if doplot>=doplot_light
                titlestring = sprintf('h-frame %d, whisker %d (%.1f)', handles.currentframe.h, w, Emin);
                set(handles.currentframe_display,'String',titlestring);
                clear titlestring
            end
            
          
            % Post-process:
            % 1. Find follicle position (intersection of polyfit to bezier with
            % folliclemask).  Do this by fitting a standard 'polyfit' polynomial to
            % the bezier points (since this extrapolates better).
            % 2. Renormalise: shift control points such that arc length from fp to
            % P0 and to P1 equals that of the initial frame.
            
%             if (handles.currentframe==handles.video.startframe) || ~handles.whisker(w).tracked(lastframe)
%                 %                 ! shouldn't this be whisker specific? 2nd priority
%                 handles.trpmtrs.s0 = [];
%                 %                 handles.trpmtrs.s0.v = [];
%             end
            
            %             [r3new, fp3, fpidx, sfp, s, sp2] = ...
            %                 postproc_v3(r3new, folliclemask, handles.trpmtrs.s0, handles.dt);
            % fp3 seems not to be used - delete?
            % ditto fpidx
            % sp2:      location along Bezier at which t=1 (ie arc-position of P2) wrt P0 (s=0)
            % sfp:      location along Bezier corresonding to fp wrt P0
            %           (fp is location of intersection between whisker Bezier and folliclemask)
            
            % 220216 for moment (debugging), ignore postproc:
            %             r3new_tmp = r3new;
            try
            [r3new, fp3, ~, sfp, s, sp2] = ...
                postproc_wm4a(r3new, folliclemask, handles.whisker(w).s0, handles.dt, handles.calib);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                    handles.whisker(w).selected=false;
                    attempts(attempt).Emin=handles.whisker(w).energy_threshold;
                    fp3=zeros(1,3);
                    r3new=zeros(3,3);
                    sfp=1;
                    sp2=1;
                    attempts(attempt).doplot = fp3;
                    attempts(attempt).r3new = r3new;
                    continue
                end
                %rethrow(ME)
            end
%             % WM4a 200216:
%             if (handles.currentframe==handles.video.startframe) || ~handles.whisker(w).tracked(lastframe)
%                 handles.trpmtrs.s0(1) = -sfp;
%                 handles.trpmtrs.s0(2) = sp2;
%             end
            % bug(?) fix 250216:
            if (handles.currentframe.h==handles.video.h.startframe) || ~handles.whisker(w).tracked(lastframe)
                handles.whisker(w).s0(1) = -sfp;
                handles.whisker(w).s0(2) = sp2;
            end            
            
            if doplot>=doplot_light
                t = [0:handles.dt:1];
                b3 = bezierval(r3new,t);
                r2new = projection2(r3new,handles.calib);
                fp2 = projection2(fp3',handles.calib);
                b2 = projection2(b3,handles.calib);
                titlestring = sprintf('h-frame %d, whisker %d (%.1f)', handles.currentframe.h, w, Emin);
                axes(handles.viewh)
                %colr=wmod
                plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod},fp2(1,1),fp2(2,1),'y.')
                plot(r2new(1,:,1),r2new(2,:,1),handles.points{wmod}),
                if ~handles.tracking2D
                axes(handles.viewv)
                plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod},fp2(1,2),fp2(2,2),'y.')
                plot(r2new(1,:,2),r2new(2,:,2),handles.points{wmod}),
                set(handles.currentframe_display,'String',titlestring);
                end
                clear titlestring
                clear b3 r2new fp2 b2 t
            end
%             r3all(handles.currentframe,:,:) = r3new;
            
            % Test if solution is good  (energy is low enough).
            %   If yes, go on to next frame
            %   If no, retrack this frame using predictor=false initial conditions.
            %   If this is the second time through, accept the least bad of the two
            %   solutions but stop tracking and return control to the user.
            
            % store current solution, in case need to retrack:
            attempts(attempt).Emin = Emin;
            attempts(attempt).fp3 = fp3;
            attempts(attempt).r3new = r3new;
            clear Emin fp3 fpidx
            
            if attempt==1
                if (attempts(1).Emin<handles.whisker(w).energy_threshold)
                    % solution good
                    break % ie out of the attempt loop
                elseif firstframe(w)
                    % no point doing repeated manual initialisation
                    break
                else
                    titlestring = sprintf('h-frame %d, whisker %d (%.1f) BAD SOLUTION, attempting rescue...',handles.currentframe.h,w,attempts(1).Emin);
                    set(handles.currentframe_display,'String',titlestring);
                    clear titlestring
                end
            end
            
        end % attempt loop
        nattempts = attempt;
        
        % Take the best (or least bad) solution amongst those attempted:
        tmp = zeros(1,nattempts);
        for a = 1:nattempts
            tmp(a) = attempts(a).Emin;
        end
        clear a
        [~,bestattempt] = min(tmp);
        clear tmp
        Emin_best = attempts(bestattempt).Emin;
        r3 = attempts(bestattempt).r3new;
        r3all(handles.currentframe.h,:,:) = r3;
        fp3 = attempts(bestattempt).fp3;
        handles.whisker(w).Emin(handles.currentframe.h) = attempts(bestattempt).Emin;
        clear bestattempt
        
        % Was the best solution good enough?
        if firstframe(w)
            % avoid repeated manual initisation
            % carry on to next frame
        elseif (Emin_best<handles.whisker(w).energy_threshold)
            % yes, carry on to next frame
        else
            % no - failed - deselect current whisker
            beep
            titlestring = sprintf('h-frame%d, whisker %d (%.1f,%.1f)',handles.currentframe.h,w,attempts(1).Emin,attempts(2).Emin);
            set(handles.currentframe_display,'String',titlestring);
            clear titlestring
            handles.whisker(w).selected = false;
            update_current_whisker_display(w, handles)
        end
        clear attempts
        % check that the whisker is in the view
        if ~w_in_view(r3, im.h, im.v, handles.dt,handles.calib)
            titlestring = sprintf('h-frame%d, whisker %d out of the view',handles.currentframe.h,w);
            set(handles.currentframe_display,'String',titlestring);
            clear titlestring
            handles.whisker(w).selected = false;
            update_current_whisker_display(w, handles)
        end  
        tstar = [0;0];  % ie compute curvature etc at P0
        %         for v = 1:2  % loop over views
        %             kappa(v) = curvature(squeeze(r(:,:,v)),tstar(v));
        %             theta(v) = base_angle(squeeze(r(:,:,v)),tstar(v));
        %         end
        clear r3
        
        % Save data in handles:
        handles.whisker(w).r3all = r3all;
        handles.whisker(w).fp3_all(handles.currentframe.h,:,:) = fp3;
        handles.whisker(w).tracked(handles.currentframe.h) = true;
        %         handles.whisker(w).kappa3_all(handles.currentframe.h,:) = kappa3;
        %         handles.whisker(w).theta3_all(handles.currentframe.h,:) = theta3;
        
        % v4b5: make kappa circle plot work for 3d
        % 261016
%         if get(handles.radiobutton_kinematics_test,'Value')
%             ttest = str2double(get(handles.edit_kinematics_test_position,'String'));
%             % t=0.5 ~ P1; t=0.0 ~ P0
%             plot_kinematics_test(handles,w,wmod,ttest)
%         end
        
        clear fp3 kappa3 theta3
        clear r3all kappa3_all theta3_all fp3_all
        
        firstframe(w) = false;
        
    end % w(hisker) loop
    
    % if tracking for all whiskers, failed, stop
    if all(~([handles.whisker(:).selected]))
        track = 0;
    end
    
    % Save results from this frame to disk:
    whisker = handles.whisker;
    trpmtrs = handles.trpmtrs;
    roi = handles.roi;
    calib = handles.calib;
    %follicle=handles.folliclemask;
    save(handles.trfname,'whisker','trpmtrs','roi','calib')
    clear trpmtrs whisker roi calib
    
    % Proceed to track next frame or stop?
    % v4b2 050916: important change to logic in following lines
%     if handles.currentframe.h ~= handles.video.h.stopframe;   %v4b1
    onwards = false;
    if track_fwd && (handles.currentframe.h~=handles.video.h.stopframe)
        % we're tracking forwards and it's not the final frame, so carry on
        onwards = true;
        handles.currentframe.h = modadd(handles.currentframe.h,1,handles.video.h.nframes);
        if ~handles.tracking2D
        handles.currentframe.v = modadd(handles.currentframe.v,1,handles.video.v.nframes);
        end
    elseif ~track_fwd && (handles.currentframe.h~=handles.video.h.startframe)
        % we're tracking backwards and it's not the first frame, so carry on
        onwards = true;
        handles.currentframe.h = modsubtract(handles.currentframe.h,1,handles.video.h.nframes);
        if ~handles.tracking2D
        handles.currentframe.v = modsubtract(handles.currentframe.v,1,handles.video.v.nframes);
        end
    end

    if onwards
    %         handles.frame = load_frame_2views(handles.video,handles.currentframe,handles.roi);%,...
        handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
        if ~handles.tracking2D
        handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);
        end
%         set(handles.currentframe_display,'string',num2str(handles.currentframe.h))
        if doplot
            pause(.01)
            cla(handles.viewh,'reset')
            titles.h = num2str(handles.currentframe.h); 
            haxes.h = handles.viewh; 
            
            if ~handles.tracking2D
            cla(handles.viewv,'reset')
            titles.v = num2str(handles.currentframe.v);
            haxes.v = handles.viewv; 
            end
            plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
            clear haxes titles
            
        end
    else
        % end of video, stop
        track = 0;
        beep,pause(.5),beep,pause(.5),beep,pause(.5),beep,pause(.5),beep
    end
    clear onwards
    
    pause(.01)
    
    %%%%%tmp debugging Andrea
    %track the first 10 frames and continue with next video
%     if handles.currentframe.h>=(handles.video.h.startframe+10)
%         current=handles.currentframe.h-handles.video.h.startframe
%         pause(0.5)
%         track=0;
%     end
    
    
end % while loop
clear track
catch
     fnoutput=2;
     disp(['Error during tracking' handles.fname.h])
     return
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function nvec = unit_normals(tvec)
% 
% n = sqrt(sum(tvec.^2,1));
% tvec = tvec ./ (ones(2,1)*n);   % unit tangent vectors
% clear n
% R = [0 -1;1 0]; % 90deg rotation matrix
% nvec = R*tvec;       % normal vectors at s=[0,1]
% n = sqrt(sum(nvec.^2,1));
% nvec = nvec ./ (ones(2,1)*n);    % unit normals
% clear tvec R n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
function in_view=w_in_view(r3, imh, imv, dt,calib)
in_view=true;

% find the points on the projected bezier curves that are internal to the images:
mv = calib.matrix(1,:);
mw = calib.matrix(2,:);
ov = calib.vector(1);
ow = calib.vector(2);

t = 0:dt:1;
b = bezierval(r3,t);
x = b(1,:);
y = b(2,:);
v = mv*b + ov;
w2 = mw*b + ow;

% find the points on the projected bezier curves that are internal to the images:
wid = size(imh.s,2);
hgt = size(imh.s,1);
w_in_H=~isempty(find((x>=1)&(x<=wid)&(y>=1)&(y<=hgt)));
%if gdtp.h

if isfield(imv,'s')
wid = size(imv.s,2);
hgt = size(imv.s,1);
w_in_V = ~isempty(find((v>=1)&(v<=wid)&(w2>=1)&(w2<=hgt)));
else
w_in_V=true;    
end
in_view=w_in_H && w_in_V;
clear wid hgt gdpt x y w2 v w_in_V w_in_H




function np = normal_plane(r,t0)
% Compute normal plane to the tangent to the Bezier curve defined by the
% control points (columns of 'r') at t==t0
%
% 'np' is 3x2: its columns span the normal plane

tv = bezierdtval(r,t0); % tangent vector
tv = tv/norm(tv);   % unit vector

% Orthogonalise the unit cartesian vectors wrt t (see notebook entry
% 151115):

M = eye(3) - tv*tv';  %This is a 3x3 matrix with rank 2

% It is the non-null subspace of M that we want:
[vec,val] = eig(M);
[~,idx] = sort(diag(val),'descend');
% extract eigenvectors with highest eigenvalue:
np = vec(:,idx(1:2));

%% windowbutton callbacks

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vpos = get(handles.viewh,'position'); % origin is upper left
pt = get(hObject,'CurrentPoint');    % mouse click location wrt figure
% set(handles.text_mouse,'String',num2str(pt))

% axes(handles.viewh)
% hold on

% location wrt viewh axes:
tpt(1) = pt(1)-vpos(1);
tpt(2) = vpos(4)-(pt(2)-vpos(2));
% plot(tpt(1),tpt(2),'bx')
% set(handles.text_axes,'String',num2str(tpt))
handles.CurrentPoint = tpt;




%% subfunctions

function [rinit, goodsolution] = auto_initialise(trfile, handles)

    if ~handles.tracking2D
        load(trfile,'-mat','whisker','calib')
        handles.calib=calib;
    else
        load(trfile,'-mat','whisker')
        %define calib matrix as dummy values in 2D version
        calib.matrix=zeros(2,3);
        calib.vector=zeros(2,1);
        calib.mv=zeros(1,3);
        calib.mw=zeros(1,3);
        calib.ov=0;
        calib.ow=0;
        handles.calib = calib;
    end
% % choose a representative subset of frames:
% [theta,frame] = sort(theta_all);
% theta = theta(2:end-1);
% frame = frame(2:end-1);
% n = length(theta);
% % Crude way to pick the subset, since biased to the most common angle.
% % But, so long as the subset is big enough, it should satisfice.
% interval = 2;%max([floor(n/100),1]);
% frame_subset = frame(1:interval:end);

% do PCA on the concatenated control point vectors:
debug =false;
if debug 
    d3=figure;
end

for w = 1:handles.Nwhiskers
    trframe{w} = find(whisker(w).tracked);
end
clear w
trframe = unique([trframe{:}]);
data = [];
for w = 1:handles.Nwhiskers
    if whisker(w).selected
        tmp = whisker(w).r3all(trframe,:,:);
        data = [data tmp(:,:)];
    end
end
clear w tmp
% [vec,val] = eig(cov(data));
% 
% if 0
%     figure
%     subplot 221
%     imagesc(cov(data))
%     subplot 222
%     plot(diag(val),'.-')
%     subplot 223
%     plot(vec(:,end:-1:end-1))
%     subplot 224
%     plot(data*vec(:,end),data*vec(:,end-1),'.')
%     pause
% end

frame.h = double(handles.frame.h(:,:,2))-handles.meanframe.h(:,:,2);
im.h.raw = frame.h;
im.h.s = imfilter(im.h.raw,handles.gaussian);
[im.h.dx,im.h.dy] = gradient(im.h.s);

if ~handles.tracking2D
    frame.v = double(handles.frame.v(:,:,2))-handles.meanframe.v(:,:,2);
    im.v.raw = frame.v;
    im.v.s = imfilter(im.v.raw,handles.gaussian);
    [im.v.dx,im.v.dy] = gradient(im.v.s);
else
     im.v=[];
end
% 
% z = data*vec(:,end);
% zrange = linspace(-2*std(z),2*std(z),100);
if debug
%     figure, 
     subplot 121, hold on
    imagesc(frame.h), colormap gray, set(gca,'ydir','reverse')
    subplot 122, hold on
    imagesc(frame.v), colormap gray, set(gca,'ydir','reverse')
end
% for zidx = 1:length(zrange)
%     z = zrange(zidx);
%     tmp = z*vec(:,end)+mean(data)';
%     Ew = zeros(1,handles.Nwhiskers);
%     for w = 1:handles.Nwhiskers
%         tmp2 = tmp(9*(w-1)+1:9*w);
%         r3 = reshape(tmp2,3,3);
%         if debug,
%             plot(r3(1,:),r3(2,:),handles.colours{w})
%         end
%         z0 = [0 0 r3(:,2)' 0 0];
%         np0 = normal_plane(r3,0);   % np0 ~ column vectors that span the plane
%         np1 = normal_plane(r3,1);
%         [Ew(w),~] = Bimfun_wm4(z0, r3, np0, np1, im.h, im.v, handles.dt, ...
%             handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);
%     end
%     E(zidx) = sum(Ew);
%     title(E(zidx)),pause
%     clear tmp tmp2 r3 z0 np0 np1 Ew    
% end
% clear zidx

testframes = round(linspace(trframe(1),trframe(end),500));
%testframes = trframe(1)

Ew = zeros(length(testframes),handles.Nwhiskers);
for fridx = 1:length(testframes)
    fr = testframes(fridx);
    for w = 1:handles.Nwhiskers
%         if ~whisker(w).selected
%             continue
%         end
        wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
        r3 = squeeze(whisker(w).r3all(fr,:,:));
        r2 = projection2(r3,handles.calib);
        if debug
%             subplot 121
%             plot(r2(1,:,1),r2(2,:,1),handles.lines{wmod})
%             subplot 122
%             plot(r2(1,:,2),r2(2,:,2),handles.lines{wmod})
        end
        z0 = [0 0 r3(:,2)' 0 0];
        np0 = normal_plane(r3,0);   % np0 ~ column vectors that span the plane
        np1 = normal_plane(r3,1);
        %[E,gradE] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib)
        [Ew(fridx,w),~] = Bimfun_wm4(z0, r3, np0, np1, im.h, im.v, handles.dt, ...
            handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);
    end
end
clear fridx fr w r3 r2 wmod

% find the contour/frame with the lowest E:
goodsolution(1:handles.Nwhiskers) = false;
EminH=zeros(1,3);
EminV=zeros(1,3);
Emin=zeros(1,3);
EwV=zeros(21,handles.Nwhiskers);
EwH=zeros(21,handles.Nwhiskers);
for w = 1:handles.Nwhiskers
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    %find the frame that minimise the set of whiskers
    %size(sum(Ew,2))
    [~,fridx] = min(sum(Ew,2));
    
    Emin(w)=Ew(fridx,w);
    %minimise every whisker independently
    %[Emin(w),fridx] = min(Ew(:,w));
    fr = testframes(fridx);
    %%best solution for whisker w
    
    rinit{w} = squeeze(whisker(w).r3all(fr,:,:));
    r2 = projection2(rinit{w},handles.calib);
    
        %%%%create 10 traslations on horizontal view
    %%%translation in horizontal view
    %%rinit{w} is the original solution
    %%rinittras{w} is the translation in x axis

    v_tras=zeros(size(rinit{w}));
    rinittras=zeros(21,3,3);
    for i=-10:10
        %%%%how much to traslate
        v_tras(1,:)=i;
        rinittras(i+11,:,:)=rinit{w}+v_tras;
    end
    %     %compute the energy for the translations
    for fridx = 1:size(rinittras,1)
        r3 = squeeze(rinittras(fridx,:,:));
        z0 = [0 0 r3(:,2)' 0 0];
        np0 = normal_plane(r3,0);  % np0 ~ column vectors that span the plane
        np1 = normal_plane(r3,1);
        %[E,gradE] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib)
        
        [EwH(fridx,w),~] = Bimfun_wm4(z0, r3, np0, np1, im.h, im.v, handles.dt, ...
            handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);    
        clear z0 np0 np1   
    end
    
    %select best translation
    [EminH(w),ntras] = min(EwH(:,w));
    %%best solution for whisker w from the second process
    rinit2{w} = squeeze(rinittras(ntras,:,:));
    r2new = projection2(rinit2{w},handles.calib);
    rinit{w}=rinit2{w};
    EminV(w)=EminH(w);
    if debug
        rinittras1{w}=squeeze(rinittras(1,:,:));
        rinittras2{w}=squeeze(rinittras(21,:,:));
        
        r2tras1=projection2(rinittras1{w},handles.calib);
        r2tras2=projection2(rinittras2{w},handles.calib);
        subplot 121
        plot(r2(1,:,1),r2(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2tras1(1,:,1),r2tras1(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2tras2(1,:,1),r2tras2(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2new(1,:,1),r2new(2,:,1),'y','linewidth',1)
        subplot 122
        plot(r2(1,:,2),r2(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2tras1(1,:,2),r2tras1(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2tras2(1,:,2),r2tras2(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2new(1,:,2),r2new(2,:,2),'y','linewidth',1)
        pause(1)
    end
   
    
    %%%%create 10 traslations on vertical view
    %%%translation in vertical view
    %%rinit{w} is the original solution
    %%rinittras{w} is the translation in z axis
    %         b3 = bezierval(rinit{w},[0:0.1:1]);
    if ~handles.tracking2D
    v_tras=zeros(size(rinit{w}));
    rinittras=zeros(21,3,3);
    for i=-10:10
        %%%%how much to traslate
        v_tras(3,:)=i;
        rinittras(i+11,:,:)=rinit{w}+v_tras;
    end
    %     %compute the energy for the translations
    for fridx = 1:size(rinittras,1)
        r3 = squeeze(rinittras(fridx,:,:));
        z0 = [0 0 r3(:,2)' 0 0];
        np0 = normal_plane(r3,0);  % np0 ~ column vectors that span the plane
        np1 = normal_plane(r3,1);
        %[E,gradE] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib)
        
        [EwV(fridx,w),~] = Bimfun_wm4(z0, r3, np0, np1, im.h, im.v, handles.dt, ...
            handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);    
        clear z0 np0 np1   
    end
    
    %select best translation
    [EminV(w),ntras] = min(EwV(:,w));
    %%best solution for whisker w from the second process
    
    rinit2{w} = squeeze(rinittras(ntras,:,:));
    r2new = projection2(rinit2{w},handles.calib);
    rinit{w}=rinit2{w};
    
    if debug
        rinittras1{w}=squeeze(rinittras(1,:,:));
        rinittras2{w}=squeeze(rinittras(21,:,:));
        
        r2tras1=projection2(rinittras1{w},handles.calib);
        r2tras2=projection2(rinittras2{w},handles.calib);
        subplot 121
        plot(r2(1,:,1),r2(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2tras1(1,:,1),r2tras1(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2tras2(1,:,1),r2tras2(2,:,1),handles.lines{wmod},'linewidth',1)
        plot(r2new(1,:,1),r2new(2,:,1),'y','linewidth',1)
        subplot 122
        plot(r2(1,:,2),r2(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2tras1(1,:,2),r2tras1(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2tras2(1,:,2),r2tras2(2,:,2),handles.lines{wmod},'linewidth',1)
        plot(r2new(1,:,2),r2new(2,:,2),'y','linewidth',1)
    end
    

    end
    disp(sprintf('Best solution for whisker %d: E=%.1f ',w,EminV(w))),
    if EminV(w)<handles.whisker(w).energy_threshold*1.1
        goodsolution(w) = true;
        disp('...good')
    else
        disp('...bad')
        beep
        % use manual initialisation
    end
end
clear w wmod fr r2 Emin EminH EminV EwV EwH


% if above didn't work, try PCA approach
% 010316: currently there is a bug in this
% if debug
%     figure, subplot 121, hold on
%     imagesc(frame.h), colormap gray, set(gca,'ydir','reverse')
%     subplot 122, hold on
%     imagesc(frame.v), colormap gray, set(gca,'ydir','reverse')
% end
% for w = 1:handles.Nwhiskers
%     wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
%     if goodsolution(w)
%         break
%     end
%     disp('Trying PCA approach...')
%     trframe = find(whisker(w).tracked);
%     data = whisker(w).r3all(trframe,:,:);
%     data = data(:,:);
%     [vec,val] = eig(cov(data));
%     z = data*vec(:,end);
%     zrange = linspace(-4*std(z),4*std(z),100);
%     Ew = zeros(1,numel(zrange));
%     for zidx = 1:length(zrange)
%         z = zrange(zidx);
%         tmp = z*vec(:,end)+mean(data)';
%         r3 = reshape(tmp,3,3);
%         r2 = projection2(r3,handles.calib);
%         if debug
%             subplot 121
%             plot(r2(1,:,1),r2(2,:,1),handles.lines{wmod})
%             subplot 122
%             plot(r2(1,:,2),r2(2,:,2),handles.lines{wmod})
%         end
%         z0 = [0 0 r3(:,2)' 0 0];
%         np0 = normal_plane(r3,0);   % np0 ~ column vectors that span the plane
%         np1 = normal_plane(r3,1);
%         [Ew(zidx),~] = Bimfun_wm4(z0, r3, np0, np1, im.h, im.v, handles.dt, ...
%             handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);
%     end
%     clear z tmp r3 r2 z0 np0 np1
%     [Emin,zidx] = min(Ew);
%     z = zrange(zidx);
%     tmp = z*vec(:,end)+mean(data)';
%     r3 = reshape(tmp,3,3);
%     clear data vec val zrange z tmp
%     rinit{w} = r3;
%     r2 = projection2(rinit{w},handles.calib);
%     if debug
%         subplot 121
%         plot(r2(1,:,1),r2(2,:,1),handles.lines{wmod},'linewidth',5)
%         subplot 122
%         plot(r2(1,:,2),r2(2,:,2),handles.lines{wmod},'linewidth',5)
%     end
%     disp(sprintf('Best solution for whisker %d: E=%.1f ',w,Emin)),
%     if Emin<handles.whisker(w).energy_threshold
%         goodsolution(w) = true;
%         disp('...good')
%     else
%         disp('...bad')
% %         beep
%     end
% end    
% clear w wmod



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function frame = load_frame(videopmtrs, framenum, roi)

switch videopmtrs.type
%     case 'avi'
%         video = read(videopmtrs.vObj, framenum);
%         frame = video(:,:,:,1);
    case 'dat'
        offset = videopmtrs.header.imagesize * (framenum-1) + videopmtrs.offset;
        fseek(videopmtrs.fid,offset,-1);
        tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
        tmp = reshape([tmp; zeros(24,1)],videopmtrs.width.raw,videopmtrs.height.raw)';
        frame = uint8(zeros(videopmtrs.height.raw,videopmtrs.width.raw,3));
        frame(:,:,1) = tmp;
        frame(:,:,2) = tmp;
        frame(:,:,3) = tmp;
        clear tmp
    otherwise
        error('Unhandled video file type')
end

frame = frame(roi(3):roi(4),roi(1):roi(2),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function frame = load_frame_2views(videopmtrs,frameidx,roi)%,rotatevview)
% 
% switch videopmtrs.type
%     case 'avi'
%         video = read(videopmtrs.vObj, frameidx);
%         frame.raw = video(:,:,:,1);
%     case 'dat'
%         offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
%         fseek(videopmtrs.fid,offset,-1);
%         tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
%         tmp = reshape([tmp; zeros(24,1)],videopmtrs.width.raw,videopmtrs.height.raw)';
%         frame.raw = uint8(zeros(videopmtrs.height.raw,videopmtrs.width.raw,3));
%         frame.raw(:,:,1) = tmp;
%         frame.raw(:,:,2) = tmp;
%         frame.raw(:,:,3) = tmp;
%         clear tmp
%     otherwise
%         error('Unhandled video file type')
% end
% 
% 
% % extract ROIs:
% frame.h = frame.raw(roi.h(3):roi.h(4),roi.h(1):roi.h(2),:);
% frame.v = frame.raw(roi.v(3):roi.v(4),roi.v(1):roi.v(2),:);
% % rotatevview:
% tmp = frame.v(:,:,1)';
% frame.v = zeros(size(tmp,1),size(tmp,2),3);
% frame.v(:,:,1) = tmp;
% frame.v(:,:,2) = tmp;
% frame.v(:,:,3) = tmp;
% clear tmp
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = predict_contour(rold, ddr, miumax)
% The predictor is a linear combination of two terms.
% Term 1: P0 & P2 come from r = rold + ddr
%         'P1prev' is set so that the contour has the same shape
%         (curvature) as in the previous frame
% Term 2: P0 & P2 as above.
%         'P1linear' is the midpoint of P0 and P2.  The point of this is
%         that, when the contour is linear, the P1 coordinate is
%         ill-defined (any point along the P0-P2 line works).
%         Using the midpoint gives more stable extrapolation
%         beyond the range t=[0,1] and avoids instability.
% P1prev and P1linear are combined according to how linear the contour
% was in the previous frame.  If it was near-linear, P1linear is
% weighted most; if it was curved, P1prev is weighted most.
% Computation of term 1:
P0old = rold(:,1);
P1old = rold(:,2);
P2old = rold(:,3);
% P0 and P2 are easy:
P0 = P0old + ddr(:,1);
P2 = P2old + ddr(:,3);
% To set P1prev, the idea is that the shape of the quadratic Bezier
% curve is defined by the triangle P0-P1-P2.  So, once we've moved
% P0 and P2, we need to move P1 such that the triangle
% P0-P1-P2 is the same size and shape (but not orientation or
% position) as the triangle P0old-P1old-P2old.
% So, compute angle of line P0-P1 wrt the line P0-P2 so that we can
% find P1prev by rotation wrt the new P0-P2 line:
tmp = ((P2old-P0old)'*(P1old-P0old))/(norm(P2old-P0old)*norm(P1old-P0old));
% dirty hack 300114 copied 04/03/2014
tmp = min([tmp 1]);
theta = acos(tmp);

n = [0 -1; 1 0]*(P2old-P0old);  %normal to P0-P2
if (P1old-P0old)'*n > 0
    % sign of angle is correct
else
    theta = -theta;
end
clear tmp n
% Set P1 by rotation wrt the (new) P0-P2 and rescaling wrt the
% length of the line P0old-P1old:
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
P1prev = P0 + norm(P1old-P0old)*R*(P2-P0)/norm(P2-P0);
clear theta R
% Computation of term 2:
P1linear = (P0+P2)/2;
% Linearly combine P1prev and P1linear.  To determine the balance
% between them, quantify how curved the curve of the previous frame
% was.
% 'miu' is length of the normal from P1 to P0-P2, divided by the
% length of P0-P2.  mui = 0 for a straight line.  miu>0 for a
% curve:
miu = norm(P0old-2*P1old+P2old)/norm(P2old-P0old);
% Introduce a sensitivity parameter.  When miu>=miumax, P1=P1prev.
P1 = min([miu/miumax,1])*P1prev + (1-min([miu/miumax,1]))*P1linear;
clear P1prev P1linear miu miumax
r = [P0 P1 P2];
clear P0old P1old P2old P0 P1 P2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = predict_contour_wm4(rold, ddr)
% Simplified predictor that doesn't use miumax
% Model: P0 & P2 come from r = rold + ddr
% see the old function ('predict_contour') for a more sophisticated approach.        
P0old = rold(:,1);
P1old = rold(:,2);
P2old = rold(:,3);
% P0 and P2 are easy:
P0 = P0old + ddr(:,1);
P2 = P2old + ddr(:,3);
% In this version (_mw4) also do basic thing for P1:
P1 = P1old + ddr(:,3);

% following commented out 250216:
% % To set P1prev, the idea is that the shape of the quadratic Bezier
% % curve is defined by the triangle P0-P1-P2.  So, once we've moved
% % P0 and P2, we need to move P1 such that the triangle
% % P0-P1-P2 is the same size and shape (but not orientation or
% % position) as the triangle P0old-P1old-P2old.
% % So, compute angle of line P0-P1 wrt the line P0-P2 so that we can
% % find P1 by rotation wrt the new P0-P2 line:
% tmp = ((P2old-P0old)'*(P1old-P0old))/(norm(P2old-P0old)*norm(P1old-P0old));
% % dirty hack 300114 copied 04/03/2014
% tmp = min([tmp 1]);
% theta = acos(tmp);
% 
% n = [0 -1; 1 0]*(P2old-P0old);  %normal to P0-P2
% if (P1old-P0old)'*n > 0
%     % sign of angle is correct
% else
%     theta = -theta;
% end
% clear tmp n
% % Set P1 by rotation wrt the (new) P0-P2 and rescaling wrt the
% % length of the line P0old-P1old:
% R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% P1 = P0 + norm(P1old-P0old)*R*(P2-P0)/norm(P2-P0);
% clear theta R

r = [P0 P1 P2];
% clear P0old P1old P2old P0 P1 P2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew, fp, fpidx, tfp, sfp, s, tstar, sp2] = postproc(r, folliclemask, s0, dt, sstar, extrap_type)
% sp2:      arc length position at which t=1 (ie arc-position of P2)
% fp:       intersection of whisker contour with folliclemask (x,y coordinate)
% fpidx:    folliclemask(fidx) is the point on folliclemask nearest to fp
% tfp:      "t" coordinate along whisker contour of fp
% sfp:      "s" coordinate of fp
% tstar:   "t" coordinate of sstar - ie the point, wrt tfp, at which to take curvature measurements

% compute arc length parameter 's':
t = -1:dt:1.5;
idx = find(t==0);
b = bezierval(r,t);
db = diff(b,1,2);
ds = sqrt(sum(db.^2,1));
s = [0 cumsum(ds)]; % s(t=-1) = 0
s = s - s(idx); % now s(t=0) = 0
idx = find(t==1);
sp2 = s(idx);
clear db ds idx

% Find intersection between whisker polynomial and folliclemask
switch extrap_type
    case 'bezier'
        [fp, fpidx, tfp] = find_follicle(r,dt,folliclemask);
    case 'poly'
        %         [fp, fpidx, tfp] = find_follicle_poly(r,dt,folliclemask);
        [fp] = find_follicle_poly(r,dt,folliclemask);
    otherwise
        error('Should never happen')
end
% if B(s) is the (x,y) point at 's', then fp = B(s(fpidx)).  this is used in snout detection.
% tfp is needed to sfp - compute 's' coordinate of the snout-whisker
% intersection, and this, in turn, is used for the contour repositioning.
[~,idx] = min((t-tfp).^2);
if length(idx)~=1
    error('')
end
sfp = s(idx);
clear idx

% Renormalise contour using polynomial (non-Bezier) approach:
% Relocate P0 and P2 so that there are at standardised distances from fp
% (s0) and adjust P1 to preserve shape.
if ~isempty(s0)
    % s0(1) is location of P0 wrt fp on the reference contour (firstframe)
    rnew = zeros(size(r));
    %             rnew(:,1) = [polyval(p(1,:),sfp+s0(1));polyval(p(2,:),sfp+s0(1))];
    %             rnew(:,3) = [polyval(p(1,:),sfp+s0(2));polyval(p(2,:),sfp+s0(2))];
    %     idx = find(t==tfp);
    sprime = s - sfp;   % now sprime(t=tfp) = 0
    d = (sprime-s0(1)).^2;
    [~,idx] = min(d);   % index of sprime ~ target P0-fp distance
    rnew(:,1) = b(:,idx);
    clear d idx
    d = (sprime-s0(2)).^2;
    [~,idx] = min(d);   % index of sprime ~ target P2-fp position
    rnew(:,3) = b(:,idx);
    clear d idx
    % find P1:
    % trying out another way 170414:
    % old way:
    %     rnew(:,2) = 2*(bezierval(r,.5) - .25*rnew(:,1) - .25*rnew(:,3));
    % new way
    rnew(:,2) = 0.5*(rnew(:,1)-r(:,1)) + 0.5*(rnew(:,3)-r(:,3)) + r(:,2);
    
else
    % it's the first frame, so using this to define s0
    rnew = r;
end

[~,idx] = min((s-sfp-sstar).^2);
tstar = t(idx);
clear idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew, fp, fpidx, sfp, s, sp2] = postproc_v3(r, folliclemask, s0, dt)
% sp2:      arc length position at which t=1 (ie arc-position of P2)
% fp:       intersection of whisker contour with folliclemask (x,y coordinate)
% fpidx:    folliclemask(fidx) is the point on folliclemask nearest to fp
% sfp:      "s" coordinate of fp
% tstar:   "t" coordinate of sstar - ie the point, wrt tfp, at which to take curvature measurements

debug = 0;

% compute arc length parameter 's':
t = -1:dt:1.5;
idx = find(t==0);
b = bezierval(r,t);
db = diff(b,1,2);
ds = sqrt(sum(db.^2,1));
s = [0 cumsum(ds)]; % s(t=-1) = 0
s = s - s(idx); % now s(t=0) = 0
idx = find(t==1);   % ie P2
sp2 = s(idx);   % length of the curve defined by 'r' between t==0 and t==1 - ie "s" location of P2 along the bezier, where s==0 ~ P0
clear db ds idx

% Find intersection between whisker polynomial and folliclemask
[fp, fpidx, p, sfp] = find_follicle_poly(r,dt,folliclemask);
% if B(s) is the (x,y) point at 's', then fp = B(s(fpidx)).  this is used in snout detection.
% tfp is needed to sfp - compute 's' coordinate of the snout-whisker
% intersection, and this, in turn, is used for the contour repositioning.
% in v3, am trying to avoid using tsp/sfp

% [~,idx] = min((t-tfp).^2);
% if length(idx)~=1
%     error('')
% end
% sfp = s(idx);
% clear idx

% Renormalise contour using polynomial (non-Bezier) approach:
% Relocate P0 and P2 so that there are at standardised distances from fp
% (s0) and adjust P1 to preserve shape.
% v3 - review the logic: P0 should be distance s0(1) along the extrapolated contour from
% fp; P2 should be s0(2) from fp.  so, compute distance along the
% extrapolated contour, starting from fp. locate P0 at the point that
% matches s0(1). located P2 at the point that matches s0(2).

if ~isempty(s0)
    % s0(1) is location of P0 wrt fp on the reference contour (firstframe)
    rnew = zeros(size(r));
    %             rnew(:,1) = [polyval(p(1,:),sfp+s0(1));polyval(p(2,:),sfp+s0(1))];
    %             rnew(:,3) = [polyval(p(1,:),sfp+s0(2));polyval(p(2,:),sfp+s0(2))];
    %     idx = find(t==tfp);
    
    % sfp is location of fp, along polynomial defined by p, where sfp==0 at P0.
    % P0:
    rnew(1,1) = polyval(p(1,:),s0(1)+sfp);
    rnew(2,1) = polyval(p(2,:),s0(1)+sfp);
    % P2:
    rnew(1,3) = polyval(p(1,:),s0(2)+s0(1)+sfp);
    rnew(2,3) = polyval(p(2,:),s0(2)+s0(1)+sfp);
    
    %     sprime = s - sfp;   % now sprime(t=tfp) = 0
    %     d = (sprime-s0(1)).^2;
    %     [~,idx] = min(d);   % index of sprime ~ target P0-fp distance
    %     rnew(:,1) = b(:,idx);
    %     clear d idx
    %     d = (sprime-s0(2)).^2;
    %     [~,idx] = min(d);   % index of sprime ~ target P2-fp position
    %     rnew(:,3) = b(:,idx);
    %     clear d idx
    % find P1:
    % trying out another way 170414:
    % old way:
    %     rnew(:,2) = 2*(bezierval(r,.5) - .25*rnew(:,1) - .25*rnew(:,3));
    % new way
    rnew(:,2) = 0.5*(rnew(:,1)-r(:,1)) + 0.5*(rnew(:,3)-r(:,3)) + r(:,2);
    
else
    % it's the first frame, so using this to define s0
    rnew = r;
end

if debug
    figure
    plot(r(1,:),r(2,:),'r.',rnew(1,:),rnew(2,:),'ro')
    %pause
    close
end
% [~,idx] = min((s-sfp-sstar).^2);
% tstar = t(idx);
% clear idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r3new, fp3, fpidx, sfp, s, sp2] = postproc_wm4a(r3, folliclemask, s0, dt, calib)
% sp2:      arc length position at which t=1 (ie arc-position of P2)
% fp:       intersection of whisker contour with folliclemask (x,y coordinate)
% fpidx:    folliclemask(fidx) is the point on folliclemask nearest to fp
% sfp:      "s" coordinate of fp

debug = 0;

% compute arc length parameter 's':
t = -1:dt:1.5;
idx = find(t==0);
b = bezierval(r3,t);
db = diff(b,1,2);
ds = sqrt(sum(db.^2,1));
s = [0 cumsum(ds)]; % s(t=-1) = 0
s = s - s(idx); % now s(t=0) = 0
idx = find(t==1);   % ie P2
sp2 = s(idx);   % length of the curve defined by 'r' between t==0 and t==1 - ie "s" location of P2 along the bezier, where s==0 ~ P0
clear db ds idx

% Find intersection between whisker polynomial and H view folliclemask
[fp3, fpidx, p, sfp] = find_follicle_poly_wm4a(r3,dt,folliclemask{1});
% if B(s) is the (x,y) point at 's', then fp = B(s(fpidx)).  this is used in snout detection.
% tfp is needed to sfp - compute 's' coordinate of the snout-whisker
% intersection, and this, in turn, is used for the contour repositioning.
% in v3, am trying to avoid using tsp/sfp

% [~,idx] = min((t-tfp).^2);
% if length(idx)~=1
%     error('')
% end
% sfp = s(idx);
% clear idx

if debug
    r2 = projection2(r3,calib);
    b2 = projection2(b,calib);
    fp2 = projection2(fp3',calib);
    figure
    subplot 221, hold on
    plot(r2(1,:,1),r2(2,:,1),'r.',fp2(1,1),fp2(2,1),'y.')
    set(gca,'ydir','reverse')
    title('h')
    subplot 222, hold on
    plot(r2(1,:,2),r2(2,:,2),'r.',fp2(1,2),fp2(2,2),'y.')
    set(gca,'ydir','reverse')
    title('v')
end
% Renormalise contour using polynomial (non-Bezier) approach:
% Relocate P0 and P2 so that there are at standardised distances from fp
% (s0) and adjust P1 to preserve shape.
% v3 - review the logic: P0 should be distance s0(1) along the extrapolated contour from
% fp; P2 should be s0(2) from fp.  so, compute distance along the
% extrapolated contour, starting from fp. locate P0 at the point that
% matches s0(1). located P2 at the point that matches s0(2).

if ~isempty(s0)
    % s0(1) is location of P0 wrt fp on the reference contour (firstframe)
    r3new = zeros(size(r3));
    %             rnew(:,1) = [polyval(p(1,:),sfp+s0(1));polyval(p(2,:),sfp+s0(1))];
    %             rnew(:,3) = [polyval(p(1,:),sfp+s0(2));polyval(p(2,:),sfp+s0(2))];
    %     idx = find(t==tfp);
    
    % sfp is location of fp, along polynomial defined by p, where sfp==0 at P0.
    % P0:
    r3new(1,1) = polyval(p(1,:),s0(1)+sfp);
    r3new(2,1) = polyval(p(2,:),s0(1)+sfp);
    r3new(3,1) = polyval(p(3,:),s0(1)+sfp);
    % P2:
    r3new(1,3) = polyval(p(1,:),s0(2)+s0(1)+sfp);
    r3new(2,3) = polyval(p(2,:),s0(2)+s0(1)+sfp);
    r3new(3,3) = polyval(p(3,:),s0(2)+s0(1)+sfp);
    
    %     sprime = s - sfp;   % now sprime(t=tfp) = 0
    %     d = (sprime-s0(1)).^2;
    %     [~,idx] = min(d);   % index of sprime ~ target P0-fp distance
    %     rnew(:,1) = b(:,idx);
    %     clear d idx
    %     d = (sprime-s0(2)).^2;
    %     [~,idx] = min(d);   % index of sprime ~ target P2-fp position
    %     rnew(:,3) = b(:,idx);
    %     clear d idx
    % find P1:
    % trying out another way 170414:
    % old way:
    %     rnew(:,2) = 2*(bezierval(r,.5) - .25*rnew(:,1) - .25*rnew(:,3));
    % new way
    r3new(:,2) = 0.5*(r3new(:,1)-r3(:,1)) + 0.5*(r3new(:,3)-r3(:,3)) + r3(:,2);
    
else
    % it's the first frame, so using this to define s0
    r3new = r3;
end

if debug
    r2 = projection2(r3new,calib);
    subplot 221
    plot(r2(1,:,1),r2(2,:,1),'ro')
    subplot 222
    plot(r2(1,:,2),r2(2,:,2),'ro')
    %pause
    close
end
% [~,idx] = min((s-sfp-sstar).^2);
% tstar = t(idx);
% clear idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = base_angle(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
%
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% theta = atan2(-dyds,dxds)*(180/pi);
% clear order pds dxds dyds

dBdt = bezierdtval(r,t);
theta = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);
clear dBdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function theta = base_angle3(r,t)
% angle of tangent vector to Bezier defined by control points {r} at point
% t
% theta(1) is in the x-y plane
% theta(2) is in the y-z plane
% theta(3) is in the x-z plane (not sure meaningful)

dBdt = bezierdtval(r,t);    % tangent vector
theta = zeros(3,length(t));
theta(3) = atan2(-dBdt(2,:),dBdt(1,:))*(180/pi);        % azimuth - whisker pointing caudal is 0'.
theta(1) = 180-atan2(-dBdt(2,:),dBdt(3,:))*(180/pi);    % elevation - whisker pointing ventral is 0'.
theta(2) = atan2(-dBdt(3,:),dBdt(1,:))*(180/pi);
clear dBdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature(r,t)

% order = size(p,2)-1;
% pds = zeros(2,order);
% pds(1,:) = p(1,1:order).*(order:-1:1);
% pds(2,:) = p(2,1:order).*(order:-1:1);
% pds2(1,:) = pds(1,1:order-1).*(order-1:-1:1);
% pds2(2,:) = pds(2,1:order-1).*(order-1:-1:1);
%
% dxds = polyval(pds(1,:),s);
% dyds = polyval(pds(2,:),s);
% d2xds2 = polyval(pds2(1,:),s);
% d2yds2 = polyval(pds2(2,:),s);
%
% kappa = dxds.*d2yds2 - dyds.*d2xds2;
% kappa = kappa ./ (dxds.^2+dyds.^2).^(3/2);
%
% clear order pds pds2 dxds dyds d2xds2 d2yds2

dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);
kappa = (dBdt(1,:).*d2Bdt2(2,:)-dBdt(2,:).*d2Bdt2(1,:)) ./ (dBdt(1,:).^2+dBdt(2,:).^2).^(3/2);
clear dBdt d2Bt2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = curvature3(r,t)
% extension of the kappa formula to 3 dimensions (wikipedia curvature page)
dBdt = bezierdtval(r,t);
d2Bdt2 = bezierdt2val(r,t);

if numel(t)>1
    error('Generalise the code!')
end

kappa = norm(cross(dBdt,d2Bdt2))/norm(dBdt).^3;

% equivalent but more cumbersome formula:
% kappa = sqrt((d2Bdt2(3,:).*dBdt(2,:)-d2Bdt2(2,:).*dBdt(3,:)).^2+(d2Bdt2(1,:).*dBdt(3,:)-d2Bdt2(3,:).*dBdt(1,:)).^2+(d2Bdt2(2,:).*dBdt(1,:)-d2Bdt2(1,:).*dBdt(2,:)).^2)./ ...
%     (dBdt(1,:).^2+dBdt(2,:).^2+dBdt(3,:).^2).^(3/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function bez = bezierfit(r,order)
% % Fit bezier curve of order 'order' to control points r
% % Initial conditions:
% switch order
%     case 1
%         bez0 = zeros(2,2);
%         bez0(:,1) = r(:,1);
%         bez0(:,2) = r(:,end);
%     case 2
%         bez0 = zeros(2,3);
%         bez0(:,1) = r(:,1);
%         bez0(:,3) = r(:,end);
%         m = ceil(size(r,2)/2);
%         bez0(:,2) = r(:,m);
%     otherwise
%         error('write more code!')
% end
% clear m
% % Minimise regression error:
% [bez,Emin,exitflag,output] = fminunc(@(bez) Bfun(bez,r),bez0,....
%                 optimset('Display','off'));%'TolX',.1,'TolFun',.01));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E,gradE] = Bimfun_old(z, r, n, im, dt, sigma)
% Assume order 2 bezier curve
bez = zeros(2,3);
bez(:,1) = r(:,1) + z(1)*n(:,1);
bez(:,2) = z(2:3)';
bez(:,3) = r(:,3) + z(4)*n(:,2);
Ereg = 0.5 * sigma * sum(sum((bez-r).^2));
gradEreg = zeros(4,1);
gradEreg(1) = sigma * (bez(:,1)-r(:,1))'*n(:,1);
gradEreg(2:3) = sigma * (bez(:,2)-r(:,2));
gradEreg(4) = sigma * (bez(:,3)-r(:,3))'*n(:,2);
clear r
t = 0:dt:1;
r = bezierval(bez,t);
% dr = diff(r,1,2);
% ds = sqrt(sum(dr.^2,1));
% s = [0 cumsum(ds)];
% snew = linspace(s(2),s(end-1),length(s));
% clear dr ds
% % figure, plot(s,t,'.',snew,0,'.'),pause
% tnew = interp1(s,t,snew);
% r = bezierval(bez,tnew);
% t = tnew;
% clear s snew tnew

% find the points on the bezier curve that are internal to the image:
w = size(im.s,2);
h = size(im.s,1);
gdpt = find((r(1,:)>=1)&(r(1,:)<=w)&(r(2,:)>=1)&(r(2,:)<=h));
clear w h

E = dt*sum(interp2(im.s,r(1,gdpt),r(2,gdpt)));

% fprintf('Eim=%.2f Ereg=%.2f Etot=%.2f\n',E,Ereg,E+Ereg)
E = E + Ereg;
gradE = zeros(4,1);
idx = interp2(im.dx,r(1,gdpt),r(2,gdpt));
idy = interp2(im.dy,r(1,gdpt),r(2,gdpt));
gradE(1) = dt*sum(idx.*((1-t(gdpt)).^2)*n(1,1)+idy.*((1-t(gdpt)).^2)*n(2,1));
gradE(2) = dt*sum(idx.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(3) = dt*sum(idy.*(2*(1-t(gdpt)).*t(gdpt)));
gradE(4) = dt*sum(idx.*(t(gdpt).^2)*n(1,2)+idy.*(t(gdpt).^2)*n(2,2));
clear gdpt idx idy
gradE = gradE + gradEreg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E,gradE] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib)
% sigma1 keeps solution close to the previous frame's solution
% sigma2 keeps P1 near the middle
%
% see notebook 190216 for explicit equations

% extension to 3d nec before debug can be used...
% debug = 0;

% column vectors that span plane normal to bezier at t==0 - r(:,1):
na = np0(:,1);
nb = np0(:,2);
% column vectors that span plane normal to bezier at t==1 - r(:,3):
nc = np1(:,1);
nd = np1(:,2);

mv = calib.matrix(1,:);
mw = calib.matrix(2,:);
ov = calib.vector(1);
ow = calib.vector(2);

% Bezier control points:
P = zeros(3,3);
P(:,1) = r(:,1) + z(1)*na + z(2)*nb;
P(:,2) = z(3:5)';
P(:,3) = r(:,3) + z(6)*nc + z(7)*nd;

% temporal contiguity constraint:
Ereg1 = 0.5 * sigma1 * sum(sum((P-r).^2));
gradEreg1 = zeros(7,1);
gradEreg1(1) = sigma1 * (P(:,1)-r(:,1))'*na;
gradEreg1(2) = sigma1 * (P(:,1)-r(:,1))'*nb;
gradEreg1(3:5) = sigma1 * (P(:,2)-r(:,2));
gradEreg1(6) = sigma1 * (P(:,3)-r(:,3))'*nc;
gradEreg1(7) = sigma1 * (P(:,3)-r(:,3))'*nd;

% Constraint to centre the middle control point:
p = P(:,3)-P(:,1);
q = P(:,2)-P(:,1);
Ereg2 = 0.5 * sigma2 * (q'*p/norm(p) - .5*norm(p)).^2;
gradEreg2 = zeros(7,1);
gradEreg2(3:5) = sigma2 * (q'*p/norm(p) - .5*norm(p)) * p;

% following will need extension to 3d...
% if debug
%     figure, subplot 121, hold on
%     plot(P(1,:),P(2,:),'r.')
%     subplot 122, hold on
%     plot(p(1,:),p(2,:),'o',q(1,:),q(2,:),'sq')
%     legend('p','q')
%     pause
%     close
% end
clear p q

clear r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the 'posterior' energy term - ie line integral of bezier over the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:dt:1;
b = bezierval(P,t);
x = b(1,:);
y = b(2,:);
v = mv*b + ov;
w = mw*b + ow;

% find the points on the projected bezier curves that are internal to the images:
wid = size(imh.s,2);
hgt = size(imh.s,1);
gdpt.h = find((x>=1)&(x<=wid)&(y>=1)&(y<=hgt));
Eh = dt*sum(interp2(imh.s,x(gdpt.h),y(gdpt.h)));

if isfield(imv,'s')
wid = size(imv.s,2);
hgt = size(imv.s,1);
gdpt.v = find((v>=1)&(v<=wid)&(w>=1)&(w<=hgt));
Ev = dt*sum(interp2(imv.s,v(gdpt.v),w(gdpt.v)));
else 
    Ev=0;
end
clear wid hgt
E = Eh + Ev + Ereg1 + Ereg2;



% compute the gradients of Eh and Ev wrt the 'z' pmtrs:
% gradients of images:
id.x = interp2(imh.dx,x(gdpt.h),y(gdpt.h));
id.y = interp2(imh.dy,x(gdpt.h),y(gdpt.h));

% d Eh/dzi:
th = t(gdpt.h);
gradEh = zeros(7,1);
gradEh(1) = dt*sum(id.x.*((1-th).^2)*na(1)+id.y.*((1-th).^2)*na(2));
gradEh(2) = dt*sum(id.x.*((1-th).^2)*nb(1)+id.y.*((1-th).^2)*nb(2));
gradEh(3) = dt*sum(id.x.*(2*(1-th).*th));
gradEh(4) = dt*sum(id.y.*(2*(1-th).*th));
gradEh(5) = 0;
gradEh(6) = dt*sum(id.x.*(th.^2)*nc(1)+id.y.*(th.^2)*nc(2));
gradEh(7) = dt*sum(id.x.*(th.^2)*nd(1)+id.y.*(th.^2)*nd(2));
clear th
% d Ev/dzi:
if isfield(imv,'dx')
id.v = interp2(imv.dx,v(gdpt.v),w(gdpt.v));
id.w = interp2(imv.dy,v(gdpt.v),w(gdpt.v));

tv = t(gdpt.v);
gradEv = zeros(7,1);
gradEv(1) = dt*sum(id.v.*((1-tv).^2)*(mv*na)+id.w.*((1-tv).^2)*(mw*na));
gradEv(2) = dt*sum(id.v.*((1-tv).^2)*(mv*nb)+id.w.*((1-tv).^2)*(mw*nb));
gradEv(3) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(1)+id.w.*(2*(1-tv).*tv)*mw(1));
gradEv(4) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(2)+id.w.*(2*(1-tv).*tv)*mw(2));
gradEv(5) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(3)+id.w.*(2*(1-tv).*tv)*mw(3));
gradEv(6) = dt*sum(id.v.*(tv.^2)*(mv*nc)+id.w.*(tv.^2)*(mw*nc));
gradEv(7) = dt*sum(id.v.*(tv.^2)*(mv*nd)+id.w.*(tv.^2)*(mw*nd));
clear tv
else
gradEv = zeros(7,1);    
end
clear gdpt id

% put it all together:
gradE = gradEh + gradEv + gradEreg1 + gradEreg2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [int] = whisker_intensity(bez, im, dt)
t = 0:dt:1;
r = bezierval(bez,t);

% find the points on the bezier curve that are internal to the image:
w = size(im.s,2);
h = size(im.s,1);
gdpt = find((r(1,:)>=1)&(r(1,:)<=w)&(r(2,:)>=1)&(r(2,:)<=h));
clear w h

int = interp2(im.s,r(1,gdpt),r(2,gdpt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = bezierval(bez,t)
% Evaluate bezier curve with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate
% function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 1
        p0 = bez(:,1);
        p1 = bez(:,2);
        b = p0 + (p1-p0)*t;
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        b = p0*(1-t).^2 + 2*p1*((1-t).*t) + p2*t.^2;
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dBdt = bezierdtval(bez,t)
% Evaluate 1st deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        %         dBdt = -p0*2*(1-t) + 2*p1*(1-2*t) + 2*p2*t;
        dBdt = 2*(p0-2*p1+p2)*t + 2*(-p0+p1)*ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dB2dt2 = bezierdt2val(bez,t)
% Evaluate 2nd deriv of bezier curve wrt t, with parameters bez(:,1),bez(:,2),... at points t
% size(t) = [1,N], where N is number of points at which to evaluate function.  Typically, t = [0,1]
% size(bez) = [2,order+1]
order = size(bez,2)-1;
switch order
    case 2
        p0 = bez(:,1);
        p1 = bez(:,2);
        p2 = bez(:,3);
        dB2dt2 = 2*(p0-2*p1+p2) * ones(size(t));
    otherwise
        error('Write more code!')
end
clear p0 p1 p2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,s] = fit_curve(r,order)
% Fit polynomial curve (parameterised by 's') to the control points:
dr = sqrt(sum(diff(r,1,2).^2,1));
s = [0 cumsum(dr)];
p = zeros(2,order+1);
p(1,:) = polyfit(s,r(1,:),order);
p(2,:) = polyfit(s,r(2,:),order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fp,fidx,tfp] = find_follicle(r,dt,folliclemask)
% fp:   intersection of whisker contour with folliclemask (snout contour)
% (x,y coordinate)
% fidx: folliclemask(fidx) is the point on folliclemask nearest to fp
% tfp:  "t" coordinate along whisker contour of the intersection point


% compute arc length parameter 's':
tneg = -.5:dt:0;
w = bezierval(r,tneg);
dw = diff(w,1,2);
ds = sqrt(sum(dw.^2,1));
% sneg = [0 cumsum(ds)] - sum(ds); % s(1)=0 ~ P0
clear dw ds

% w = zeros(2,length(tneg));
% w(1,:) = polyval(p(1,:),sneg);
% w(2,:) = polyval(p(2,:),sneg);

dist = zeros(length(tneg),size(folliclemask,2));
for i = 1:length(tneg)
    for j = 1:size(folliclemask,2)
        dist(i,j) = sum((w(:,i)-folliclemask(:,j)).^2,1);
    end
end
clear i j
[m,fidx] = min(dist');  % min over folliclemask for fixed s
% locate first minimum (in order of decreasing t)
% sidx = find(diff(m)>0,1);
[~,sidx] = min(m);      % min of the min (over s)
%     fidx = fidx(sidx);      % pt on f closest to w
% now we have Order(0) follicle position estimate:
fp0 = [w(1,sidx); w(2,sidx)];
clear dist m tmp

debug = 1;

if debug
    % plot stuff
    figure, hold on
    plot(folliclemask(1,:),folliclemask(2,:),'g')
    plot(w(1,:),w(2,:),'r',r(1,:),r(2,:),'r.')
    plot(fp0(1),fp0(2),'*')
end

% use linear interpolation to increase accuracy of fp estimate.
% find line that locally fits folliclemask (p2) and line that locally fits
% whisker contour (p1):
sidx = min([size(w,2)-1 sidx]);
if sidx>1
    p1 = polyfit(w(1,sidx-1:sidx+1),w(2,sidx-1:sidx+1),1);
    fidx = fidx(sidx);
    p2 = polyfit(folliclemask(1,fidx-1:fidx+1),folliclemask(2,fidx-1:fidx+1),1);
    % now solve for the point where p1 and p2 intersect:
    M = [p1(1) -1; p2(1) -1];
    fp = -inv(M)*[p1(2);p2(2)];
else
    fp = fp0;
end
if debug
    plot(w(1,sidx-1:sidx+1),polyval(p1,w(1,sidx-1:sidx+1)),'m.-')
    plot(folliclemask(1,fidx-1:fidx+1),polyval(p2,folliclemask(1,fidx-1:fidx+1)),'b.-')
    plot(fp(1),fp(2),'sq')
end

% deal with cases where the above interpolation fails (eg when line is
% vertical):
if (sqrt(sum(fp0-fp).^2) > 2) || sum(isinf(fp))
    % sign that the interpolation has failed
    % use the initial estimate
    fp = fp0;
end
% comment: wouldn't it be better to do the interpolation using a more
% robust, eg bezier, line representation.  ie not to make the assumption
% that y is a function of x.

clear p1 p2 M

% find the 't' value of fp.  Do this by minimising ||f(t)-fp||^2 wrt t,
% starting from initial condition 't0':
t0 = tneg(sidx);
options = optimset('display','off');
tfp = fminunc(@(t) tmpfun(t, r, fp), t0, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fp, fidx, p, sfp] = find_follicle_poly(r,dt,folliclemask)
% function [fp,fidx] = find_follicle_poly(r,dt,folliclemask)
% fp:   intersection of whisker contour with folliclemask (snout contour)
% (x,y coordinate)
% use polynomial extrapolation, not bezier ('_poly') version
% fidx: folliclemask(fidx) is the point on folliclemask nearest to fp
% p:    size(2,order+1). best-fitting polynomial approximation to the bezier
% sfp:  location, along the polynomial defined by p, of fp

debug = 0;
order = 2;

% compute arc length parameter 's':
tpos = 0:dt:1;
tneg = -1:dt:0;
wpos = bezierval(r,tpos);
wneg = bezierval(r,tneg);
dw = diff(wpos,1,2);
ds = sqrt(sum(dw.^2,1));
spos = [0 cumsum(ds)]; % from P0 to P2.  spos==0 ~ P0
dw = diff(wneg,1,2);
ds = sqrt(sum(dw.^2,1));
sneg = [0 cumsum(ds)] - sum(ds); % from t==-1 to P0.  spos==0 ~ P0
% s0idx = ceil(length(s)/2);  % assumes range of t is symmetrical about 0
clear dw ds
% s = [sneg spos(2:end)];

p(1,:) = polyfit(spos,wpos(1,:),order);
p(2,:) = polyfit(spos,wpos(2,:),order);

if debug
    % check so for
    figure, hold on
    plot(wpos(1,:),wpos(2,:),'r',r(1,:),r(2,:),'r.')
    plot(wneg(1,:),wneg(2,:),'r')
    plot(polyval(p(1,:),spos),polyval(p(2,:),spos),'m--')
    plot(polyval(p(1,:),sneg),polyval(p(2,:),sneg),'b--')
    set(gca,'ydir','reverse')
    %     pause
end

dist = zeros(length(sneg),size(folliclemask,2));
for i = 1:length(sneg)
    for j = 1:size(folliclemask,2)
        x = polyval(p(1,:),sneg(i));
        y = polyval(p(2,:),sneg(i));
        dist(i,j) = sum(([x;y]-folliclemask(:,j)).^2,1);
    end
end
clear i j x y
[m,fidx] = min(dist');  % min over folliclemask for fixed s
% locate first minimum (in order of decreasing t)
% sidx = find(diff(m)>0,1);
[~,sidx] = min(m);      % min of the min (over s)
%     fidx = fidx(sidx);      % pt on f closest to w
% now we have Order(0) follicle position estimate:
% fp0 = [w(1,sidx); w(2,sidx)];
fp0 = [polyval(p(1,:),sneg(sidx)),polyval(p(2,:),sneg(sidx))];
clear dist m tmp

if debug
    % plot stuff
    %     figure, hold on
    plot(folliclemask(1,:),folliclemask(2,:),'g')
    %     plot(w(1,:),w(2,:),'r',r(1,:),r(2,:),'r.')
    plot(fp0(1),fp0(2),'*')
    %     pause
end

% use linear interpolation to increase accuracy of fp estimate.
% find line that locally fits folliclemask (p2) and line that locally fits
% whisker contour (p1):
% p1 is tangent to the polynomial approximating the whisker at point fp0
p1(1,2) = polyval(p(1,:),sneg(sidx));
p1(1,1) = 2*p(1,1) + p(1,2);
p1(2,2) = polyval(p(2,:),sneg(sidx));
p1(2,1) = 2*p(2,1) + p(2,2);

fidx = fidx(sidx);
p2(1,1) = folliclemask(1,fidx+1)-folliclemask(1,fidx-1);
p2(1,2) = folliclemask(1,fidx-1);
p2(2,1) = folliclemask(2,fidx+1)-folliclemask(2,fidx-1);
p2(2,2) = folliclemask(2,fidx-1);

% find intersection of the p1 & p2 lines
a = p1(:,2)-p2(:,2);
M = [-p1(:,1) p2(:,1)];
s0 = inv(M)*a;
fp(1) = polyval(p1(1,:),s0(1));
fp(2) = polyval(p1(2,:),s0(1));
clear a M b

% sfp is computed to be location of fp wrt P0 along the polynomial defined
% by p; such that sfp==0 ~ P0 and sfp>0 ~ movement into the snout
% s0(1) is location of fp wrt fp0 along pdbp
% sneg(sidx) is location of fp0 along pdbp, where s=0 is at P0
sfp = s0(1) + sneg(sidx);
% sfp = - sneg(sidx);

if debug
    ds = 0:-1:-10;
    plot(polyval(p1(1,:),ds),polyval(p1(2,:),ds),'c--')
    ss = 0:.01:1;
    plot(polyval(p2(1,:),ss),polyval(p2(2,:),ss),'g--')
    plot(fp(1),fp(2),'k*')
    plot(polyval(p(1,:),sfp),polyval(p(2,:),sfp),'ksq')
    %pause
    close
end

% sidx = min([size(w,2)-1 sidx]);
% if sidx>1
% %     p1 = polyfit(w(1,sidx-1:sidx+1),w(2,sidx-1:sidx+1),1);
% %     p2 = polyfit(folliclemask(1,fidx-1:fidx+1),folliclemask(2,fidx-1:fidx+1),1);
%     % now solve for the point where p1 and p2 intersect:
%     M = [p1(1) -1; p2(1) -1];
%     fp = -inv(M)*[p1(2);p2(2)];
% else
%     fp = fp0;
% end
% if debug
%     plot(w(1,sidx-1:sidx+1),polyval(p1,w(1,sidx-1:sidx+1)),'m.-')
%     plot(folliclemask(1,fidx-1:fidx+1),polyval(p2,folliclemask(1,fidx-1:fidx+1)),'b.-')
% %     plot(fp(1),fp(2),'sq')
% end

% % deal with cases where the above interpolation fails (eg when line is
% % vertical):
% if (sqrt(sum(fp0-fp).^2) > 2) || sum(isinf(fp))
%     % sign that the interpolation has failed
%     % use the initial estimate
%     fp = fp0;
% end

% % find the 't' value of fp.  Do this by minimising ||f(t)-fp||^2 wrt t,
% % starting from initial condition 't0':
% t0 = tneg(sidx);
% options = optimset('display','off');
% tfp = fminunc(@(t) tmpfun(t, r, fp), t0, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fp3, fidx, p, sfp] = find_follicle_poly_wm4a(r3,dt,folliclemask)
% function [fp,fidx] = find_follicle_poly(r,dt,folliclemask)
% fp:   intersection of whisker contour with folliclemask (snout contour)
% (x,y coordinate)
% use polynomial extrapolation, not bezier ('_poly') version
% fidx: folliclemask(fidx) is the point on folliclemask nearest to fp
% p:    size(2,order+1). best-fitting polynomial approximation to the bezier
% sfp:  location, along the polynomial defined by p, of fp

debug = 0;
order = 2;

% compute arc length parameter 's':
tpos = 0:dt:1;
tneg = -1.5:dt:0;
wpos = bezierval(r3,tpos);
wneg = bezierval(r3,tneg);
dw = diff(wpos,1,2);
ds = sqrt(sum(dw.^2,1));
spos = [0 cumsum(ds)]; % from P0 to P2.  spos==0 ~ P0
dw = diff(wneg,1,2);
ds = sqrt(sum(dw.^2,1));
sneg = [0 cumsum(ds)] - sum(ds); % from t==-1 to P0.  spos==0 ~ P0
% s0idx = ceil(length(s)/2);  % assumes range of t is symmetrical about 0
clear dw ds
% s = [sneg spos(2:end)];

p(1,:) = polyfit(spos,wpos(1,:),order);
p(2,:) = polyfit(spos,wpos(2,:),order);
p(3,:) = polyfit(spos,wpos(3,:),order);    % new in wm4a

if debug
    % check so for
    figure, hold on
    plot(wpos(1,:),wpos(2,:),'r',r3(1,:),r3(2,:),'r.')
    plot(wneg(1,:),wneg(2,:),'r')
    plot(polyval(p(1,:),spos),polyval(p(2,:),spos),'m--')
    plot(polyval(p(1,:),sneg),polyval(p(2,:),sneg),'b--')
    set(gca,'ydir','reverse')
    %     pause
end

dist = zeros(length(sneg),size(folliclemask,2));
for i = 1:length(sneg)
    for j = 1:size(folliclemask,2)
        x = polyval(p(1,:),sneg(i));
        y = polyval(p(2,:),sneg(i));
        dist(i,j) = sum(([x;y]-folliclemask(:,j)).^2,1);    %distance in the H projection
    end
end
clear i j x y
[m,fidx] = min(dist');  % min over folliclemask for fixed s
% locate first minimum (in order of decreasing t)
% sidx = find(diff(m)>0,1);
[~,sidx] = min(m);      % min over s.
% sneg(sidx) is location along bezier closest to folliclemask - ie intersection
%     fidx = fidx(sidx);      % pt on f closest to w
% now we have Order(0) follicle position estimate:
% fp0 = [w(1,sidx); w(2,sidx)];
s0_0 = sneg(sidx);
fp0 = [polyval(p(1,:),s0_0),polyval(p(2,:),s0_0),polyval(p(3,:),s0_0)];  %x,y location of intersection in H projection
clear dist m tmp

if debug
    % plot stuff
    %     figure, hold on
    plot(folliclemask(1,:),folliclemask(2,:),'g')
    %     plot(w(1,:),w(2,:),'r',r(1,:),r(2,:),'r.')
    plot(fp0(1),fp0(2),'b*')
    %     pause
end

% use linear interpolation to increase accuracy of fp estimate.
% find line that locally fits folliclemask (p2) and line that locally fits
% whisker contour (p1):
% p1 is tangent to the polynomial approximating the whisker at point fp0
p1(1,2) = polyval(p(1,:),s0_0);
p1(1,1) = 2*p(1,1) + p(1,2);
p1(2,2) = polyval(p(2,:),s0_0);
p1(2,1) = 2*p(2,1) + p(2,2);
p1(3,2) = polyval(p(3,:),s0_0);    % new in wm4a
p1(3,1) = 2*p(3,1) + p(3,2);

fidx = fidx(sidx);
p2(1,1) = folliclemask(1,fidx+1)-folliclemask(1,fidx-1);
p2(1,2) = folliclemask(1,fidx-1);
p2(2,1) = folliclemask(2,fidx+1)-folliclemask(2,fidx-1);
p2(2,2) = folliclemask(2,fidx-1);

% find intersection of the p1 & p2 lines (still in the H projection)
a = p1(1:2,2)-p2(1:2,2);
M = [-p1(1:2,1) p2(1:2,1)];
s0 = inv(M)*a;
clear a M b

fp3(1) = polyval(p1(1,:),s0(1));
fp3(2) = polyval(p1(2,:),s0(1));
fp3(3) = polyval(p1(3,:),s0(1));    % new in wm4a

% sfp is computed to be location of fp wrt P0 along the polynomial defined
% by p; such that sfp==0 ~ P0 and sfp>0 ~ movement into the snout
% s0(1) is location of fp wrt fp0 along pdbp
% sneg(sidx) is location of fp0 along pdbp, where s=0 is at P0
sfp = s0(1) + s0_0;
% sfp = - sneg(sidx);

% the interpolation method can fall over if it happens that the two lines
% are near parallel.  if the interpolated solution is distant to the basic
% solution, use the latter:
d = norm(fp3-fp0); 
if d > 3    % units of pixels
    fp3 = fp0;
    sfp = s0_0;
end

if debug
    ds = 0:-1:-10;
    plot(polyval(p1(1,:),ds),polyval(p1(2,:),ds),'c--')
    ss = 0:.01:1;
    plot(polyval(p2(1,:),ss),polyval(p2(2,:),ss),'g--')
    plot(fp3(1),fp3(2),'k*')
    plot(polyval(p(1,:),sfp),polyval(p(2,:),sfp),'ksq')
    title(sprintf('distance between solutions %.1f',d))
    %pause
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = tmpfun(t,r,fp)
fphat = zeros(2,1);
% fphat(1) = polyval(p(1,:),s);
% fphat(2) = polyval(p(2,:),s);
fphat = bezierval(r,t);
E = sum((fp-fphat).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [folliclemask] = snout_segment(image, roi, folliclemask_old, fpidx_old, theta, snout_orientation, handles)
% snout_orientation = 'horizontal'/'vertical'

% Filter and differentiate:

ims.raw = image(roi.y,roi.x,1);
ims.w = size(ims.raw,2);
ims.h = size(ims.raw,1);
ims.mfilt = medfilt2(ims.raw,[5 5]);
% sigma = 6;
sigma = str2double(get(handles.edit_snout_sigma,'string'));
hsize = 3*[sigma sigma];
filter = fspecial('gaussian', hsize, sigma);
ims.gfilt = imfilter(ims.mfilt,filter,'replicate');
clear sigma hsize filter
% ims.Del1.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,1,2)];
% ims.Del1.y = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,1,1)];
% ims.Del2.x = [zeros(size(ims.gfilt,1),1) diff(ims.gfilt,2,2) zeros(size(ims.gfilt,1),1)];
% ims.Del2.x = [zeros(1,size(ims.gfilt,2)); diff(ims.gfilt,2,1); zeros(1,size(ims.gfilt,2))];
[ims.Del1.x,ims.Del1.y] = gradient(ims.gfilt);

% Specify the desired gradient measure:
if ~isempty(folliclemask_old)
    % Revised code (290414)
    % Go back to using normal to snout mask instead of tangent to whisker
    % solution:
    %     % use normal to snout mask in previous frame to define gradient
    %     % direction:
    tvec = folliclemask_old(:,fpidx_old) - folliclemask_old(:,fpidx_old+1);
    tvec = tvec / norm(tvec);
    R = [0 1;-1 0];
    nvec = R*tvec;
    clear tvec R
    % So following lines (290414) commented out:
    %     % use tangent to whisker solution in previous frame to define gradient
    %     % direction
    %     nvec = [-cosd(theta); sind(theta)];
    %     nvec = nvec/norm(nvec);
    
else
    % no previous frame, so use default:
    switch snout_orientation
        case 'horizontal'
            nvec = [1; 1]/sqrt(2);
        case 'vertical'
            nvec = [1; 0]/sqrt(2);
    end
end
ims.grad = nvec(1)*ims.Del1.x+nvec(2)*ims.Del1.y;

switch snout_orientation
    case 'horizontal'
        snout.x = 1:length(roi.x);
        %Andrea should fix this
        ims.grad(1:270,:)=0.1;
        [~,snout.y] = min(ims.grad);
    case 'vertical'
         ims.grad(:,1:150)=-0.2;
        snout.y = 1:length(roi.y);
        [~,snout.x] = min(ims.grad,[],2);
    otherwise
        error('Should not happen')
end

% Take every Nth pixel:
N = 5;
snout.x = snout.x(1:N:end);
snout.y = snout.y(1:N:end);
Npts = length(snout.x);
clear N

%Quality control on snout detection, debugging:
% figure
% subplot(2,2,1), imagesc(ims.raw), colormap gray
% hold on, plot(snout.x,snout.y,'r.-')
% subplot(2,2,2), imagesc(ims.gfilt)
% subplot(2,2,3), imagesc(ims.mfilt)
% subplot(2,2,4), imagesc(ims.grad)
% pause(1)

% The above produces a snout estimate that is outside the face
% slightly.  Warp the contour in a bit in direction of local normals

z = 5;
folliclemask(1,:) = snout.x + z*nvec(1,:);
folliclemask(2,:) = snout.y + z*nvec(2,:);
clear nvec

%  Debug:
%     figure
%     imagesc(ims.raw), colormap gray, colorbar
%     hold on,
%     plot(folliclemask(1,:),folliclemask(2,:),'m-',trackmask(1,:),trackmask(2,:),'g-')

% Align masks to complete image:
folliclemask(1,:) = folliclemask(1,:) + roi.x(1);
folliclemask(2,:) = folliclemask(2,:) + roi.y(1);
clear roi

% new in v4b4:
% Follicle-finding (postproc function) can fail if the folliclemask has
% outlier edges, since it can happen that when the bezier is extrapolated,
% it intersects with the outlier part of the folliclemask more closely than
% it does the "true" bit of folliclemask (see the videos noted under v4b2
% comments in the txt file and set debug==1 in the find_follicle function)
% Attempt to fix this by allowing user to set the number of extreme poitns
% to cut...

L = str2double(get(handles.edit_snout_outliers,'String'));
folliclemask(:,1:L) = [];
folliclemask(:,end-L+1:end) = [];
clear L


        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

framenum = handles.currentframe;
frameidx = framenum.h;
clear framenum

change_frame = 0;
change_cpt = 0;
change_view = 0;
change_whisker = 0;
nudge = 0;
dnudge = 0.5;
cview = handles.current_view;
supernudge = false;

switch numel(eventdata.Modifier)
    case 0
%       No modifiers:
%           horizontal arrows change frame
%           vertical arrows change control point
%           w changes whisker
        switch eventdata.Key
            case 'leftarrow'
                frameidx = modsubtract(frameidx,1,handles.video.h.nframes);
                change_frame = 1;
            case 'rightarrow'
                frameidx = modadd(frameidx,1,handles.video.h.nframes);
                change_frame = 1;
            case 'uparrow'
                cpt = handles.current_pt;
                cpt = rem(cpt+1,handles.n_control_pts+1);
                change_cpt = 1;
            case 'downarrow'
                cpt = handles.current_pt;
                cpt = (cpt-1) + (cpt==0)*(handles.n_control_pts+1);
                change_cpt = 1;
            case 'w'
                handles.current_whisker = modadd(handles.current_whisker,1,handles.Nwhiskers);
                change_whisker = 1;
        end
    case 1
%       One modifier:
%           shift/control alter rate of changing frame
%           alt enables nudging of current control point
        switch eventdata.Modifier{1}
            case 'shift'
                switch eventdata.Key
                    case 'leftarrow'
                        frameidx = modsubtract(frameidx,10,handles.video.h.nframes);
                        change_frame = 1;
                    case 'rightarrow'
                        frameidx = modadd(frameidx,10,handles.video.h.nframes);
                        change_frame = 1;
                end
            case 'control'
                switch eventdata.Key
                    case 'leftarrow'
                        frameidx = modsubtract(frameidx,100,handles.video.h.nframes);
                        change_frame = 1;
                    case 'rightarrow'
                        frameidx = modadd(frameidx,100,handles.video.h.nframes);
                        change_frame = 1;
        end
            case 'alt'
                switch eventdata.Key
                    % left/right changes x coord of CP (horizontal axis in
                    % H view - rostro-caudal)
                    case 'leftarrow'
                        nudge = -1;
                    case 'rightarrow'
                        nudge = 1;
                    % up/down changes y coord of CP (vertical axis in H
                    % view, horizonal axis in V view - medio-lateral)
                    case 'uparrow'
                        nudge = -2;
                    case 'downarrow'
                        nudge = 2;
                    % pageup/down changes z coord of CP (vertical axis in
                    % V view - dorso-ventral)
                    case 'pageup'
                        nudge = 3;
                    case 'pagedown'
                        nudge = -3;
                end
        end
    case 2
        % 'supernudge' means searching for best contour fit over small
        % range in the chosen direction
        if strmatch(eventdata.Modifier{1},'control') && strmatch(eventdata.Modifier{2},'alt')
            supernudge = true;
            switch eventdata.Key
                case 'leftarrow'
                    nudge = -1;
                case 'rightarrow'
                    nudge = 1;
                case 'uparrow'
                    nudge = -2;
                case 'downarrow'
                    nudge = 2;
                case 'pageup'
                    nudge = 3;
                case 'pagedown'
                    nudge = -3;
            end
        end
end

framenum.h = frameidx;
haxes.h = handles.viewh; 
handles.current_view = cview;

if ~handles.tracking2D
framenum.v = handles.framenums_vfromh(frameidx);
haxes.v = handles.viewv; 
end

if change_frame
    % plot new image in the baseframe.  if it exists, superimpose tracking solution.
    %     display_frame(frameidx)
%     handles.frame = load_frame_2views(handles.video,frameidx,handles.roi);%,get(handles.checkbox_rotatevview,'value'));
    titles.h = sprintf('Frame %d',framenum.h);
    if ~handles.tracking2D
    titles.v = sprintf('Frame %d',framenum.v);
    end
    if isfield(handles.whisker,'tracked')
        if handles.whisker(handles.current_whisker).tracked(frameidx) && isfield(handles.whisker,'Emin')
            titles.h = [titles.h sprintf('(%.1f)',handles.whisker(handles.current_whisker).Emin(frameidx))];
        end
    end
    handles.frame = load_and_plot_frames(handles.video,framenum,handles.roi,handles.meanframe,titles,haxes);  
    handles.currentframe = framenum;

    if isfield(handles.whisker,'r3all')
        plot_contours(framenum,handles)
    end    

        if get(handles.radiobutton_plot_snout_contour,'value') && isfield(handles,'folliclemask')&& size(handles.folliclemask,1)>frameidx
            frame_is_tracked = ~isempty(handles.folliclemask{frameidx,1});

            if frame_is_tracked
                axes(haxes.h)
                plot(handles.folliclemask{frameidx,1}(1,:),handles.folliclemask{frameidx,1}(2,:),'y-')
                if ~handles.tracking2D
                axes(haxes.v)
                plot(handles.folliclemask{frameidx,2}(1,:),handles.folliclemask{frameidx,2}(2,:),'y-')
                end
           end
        end
        
        
%%%%plot for sagital plane deleted (Andrea)
%%%%%%%%


%         for w = 1:handles.Nwhiskers
%             wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
%             if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx)
%                 %         hold on
%                 %         rall = handles.rall.h;
%                 dt = 0.05;
%                 t = 0:dt:1;
%                 r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
%                 b3 = bezierval(r3,t);
%                 r2 = projection2(r3,handles.calib);
%                 b2 = projection2(b3,handles.calib);
%                 fp3 = handles.whisker(w).fp3_all(frameidx,:)';
%                 fp2 = projection2(fp3,handles.calib);
%                 axes(handles.viewh)
%                 plot(b2(1,:,1),b2(2,:,1),handles.lines{w},r2(1,:,1),r2(2,:,1),handles.points{w},...
%                     fp2(1,1),fp2(2,1),'y.')
%                 axes(handles.viewv)
%                 plot(b2(1,:,2),b2(2,:,2),handles.lines{w},r2(1,:,2),r2(2,:,2),handles.points{w},...
%                     fp2(1,2),fp2(2,2),'y.')
%             end
%         end
        
%         keyboard
%         if get(handles.radiobutton_kinematics_test,'Value')
%             ttest = str2double(get(handles.edit_kinematics_test_position,'String'));
%             for w = 1:handles.Nwhiskers
%                 wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
%                 plot_kinematics_test(handles,w,wmod,ttest)
%             end
%             clear ttest w w mod
%         end
        
        % how do manual editing of control points?
        % could use up/downarrows to cycle between control points.  default is no
        % selection.  return key -> ginput -> updated control point. update display
%         for w = 1:handles.Nwhiskers
%             wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
%             if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx)
%                 r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
%                 axes(handles.axes_misc_plot) %#ok<LAXES>
% %                 theta = base_angle3(r3,0);
%                 tv = bezierdtval(r3,0);
%                 tv = tv/norm(tv);
%                 cv = bezierdt2val(r3,0);
%                 xpr = (eye(3)-tv*tv')*[1 0 0]';
%                 zpr = (eye(3)-tv*tv')*[0 0 1]';
%                 curv(1) = cv'*xpr;
%                 curv(2) = cv'*zpr;
%                 %                 compass(curv(1),curv(2),handles.lines{wmod})
%                 %                 axis(100*[-1 1 -1 1])
%                 compass(tv(1),tv(3),handles.lines{wmod})
%                 axis([-1 1 -1 1])
%                 title('projection of unit tangent vector at base in sagittal plane')
%                 clear dt t r3 b3 r2 b2 theta tv cv xpr zpr curv
%             end
%         end
%         clear w
        
    
    
    
end

if change_cpt
    r3 = squeeze(handles.whisker(handles.current_whisker).r3all(frameidx,:,:));
    dt = 0.02;
    t = 0:dt:1;
    b3 = bezierval(r3,t);
    %     r2 = zeros(2,3,2);
    %     r2(:,:,1) = r3(1:2,:);
    %     r2(:,:,2) = handles.calib.matrix*r3 + handles.calib.vector*ones(1,3);
    %     b2 = zeros(2,size(b3,2),2);
    %     b2(:,:,1) = b3(1:2,:);
    %     b2(:,:,2) = handles.calib.matrix*b3 + handles.calib.vector*ones(1,size(b3,2));
    r2 = projection2(r3,handles.calib);
    b2 = projection2(b3,handles.calib);
    
    %     switch cview
    %         case 1
    axes(handles.viewh), cla
    frame = double(handles.frame.h(:,:,1))-handles.meanframe.h(:,:,1);
    imagesc(frame),
    colormap gray
    hold on
    plot(b2(1,:,1),b2(2,:,1),'r-',r2(1,:,1),r2(2,:,1),'r.')
    if cpt~=0
        plot(r2(1,cpt,1),r2(2,cpt,1),'ow')
    end
    %         case 2
    if ~handles.tracking2D
    axes(handles.viewv), cla
    frame = double(handles.frame.v(:,:,1))-handles.meanframe.v(:,:,1);
    imagesc(frame),
    colormap gray
    hold on
    plot(b2(1,:,2),b2(2,:,2),'r-',r2(1,:,2),r2(2,:,2),'r.')
    if cpt~=0
        plot(r2(1,cpt,2),r2(2,cpt,2),'ow')
    end
    end
    %             r = squeeze(handles.whisker(handles.current_whisker).r3all(handles.frameidx,:,:));
    %         otherwise
    %             error('This should not happen')
    %     end
    %     size(handles.rall),handles.frameidx,pause
    %     r = squeeze(handles.rall(handles.frameidx,:,:));
    handles.current_pt  = cpt;
    clear r t b
end

% if change_view
%     switch cview
%         case 1
%             axes(handles.viewh)
%             title(sprintf('%d',handles.frameidx),'color',[1 0 0])
%             axes(handles.viewv)
%             title(sprintf('%d',handles.frameidx),'color','k')
%         case 2
%             axes(handles.viewv)
%             title(sprintf('%d',handles.frameidx),'color',[1 0 0])
%             axes(handles.viewh)
%             title(sprintf('%d',handles.frameidx),'color','k')
%         otherwise
%             error('This should not happen')
%     end
%
% end

if change_whisker
    w = handles.current_whisker;
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    titles.h = ''; titles.v = '';
    plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
    if handles.whisker(w).tracked(frameidx)
        dt = 0.05;
        t = 0:dt:1;
        r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
        b3 = bezierval(r3,t);
        r2 = projection2(r3,handles.calib);
        b2 = projection2(b3,handles.calib);
        fp2 = projection2(handles.whisker(w).fp3_all(frameidx,:)',handles.calib);
        axes(handles.viewh)
        plot(b2(1,:,1),b2(2,:,1),handles.lines{w},r2(1,:,1),r2(2,:,1),handles.points{w},...
            fp2(1,1),fp2(2,1),'y.')
        if ~handles.tracking2D
        axes(handles.viewv)
        plot(b2(1,:,2),b2(2,:,2),handles.lines{w},r2(1,:,2),r2(2,:,2),handles.points{w},...
            fp2(1,2),fp2(2,2),'y.')
        end
    end
    update_current_whisker_display(w, handles)
    titlestring = ['Whisker: ' handles.whisker(w).label ' Frame= ' num2str(frameidx)];
    if isfield(handles.whisker,'tracked')
        if handles.whisker(w).tracked(frameidx) && isfield(handles.whisker,'Emin')
            titlestring = ['Whisker: ' handles.whisker(w).label ' Frame= ' num2str(frameidx) sprintf(' Emin= (%.1f)',handles.whisker(w).Emin(frameidx))];
        end

    end
    set(handles.currentframe_display,'string',titlestring)
    clear titlestring
    clear w
end

if nudge ~= 0
    cpt = handles.current_pt;
    r3 = squeeze(handles.whisker(handles.current_whisker).r3all(frameidx,:,:,1));
%     r2 = zeros(2,3,2);
%     r2(:,:,1) = r3(1:2,:);
%     r2(:,:,2) = handles.calib.matrix*r3 + handles.calib.vector*ones(1,3);
    r2 = projection2(r3,handles.calib);
    %     sigma_prior = handles.trpmtrs.sigma_prior;
    %     sigma_prior2 = handles.trpmtrs.sigma2_prior;
    %     snout_orientation = 'horizontal';
    s0 = handles.whisker(handles.current_whisker).s0;
    %     switch cview
    %         case 1
    
    axes(handles.viewh)
    frame.h = double(handles.frame.h(:,:,1))-handles.meanframe.h(:,:,1);
    cla
    imagesc(frame.h),
    colormap gray
%     title(frameidx)
    hold on
    plot(r2(1,:,1),r2(2,:,1),'r.',r2(1,cpt,1),r2(2,cpt,1),'wo')
    %             sstar = handles.trpmtrs.sstar.h;
    %         case 2
    if ~handles.tracking2D
    axes(handles.viewv)
    frame.v = double(handles.frame.v(:,:,2))-handles.meanframe.v(:,:,2);
    cla
    imagesc(frame.v),
    colormap gray
%     title(frameidx)
    hold on
    plot(r2(1,:,2),r2(2,:,2),'r.',r2(1,cpt,2),r2(2,cpt,2),'wo')
    %             sstar = handles.trpmtrs.sstar.v;
    %         otherwise
    %             error('This should not happen')
    %     end
    %     order = handles.pmtrs.order;
    end
    r3new = r3;
    
    if supernudge
        dnudges = dnudge * (-20:20);
    else
        dnudges = dnudge;
    end
    
    r3test = r3(:,:,ones(1,numel(dnudges)));
    switch abs(nudge)
        case 1
            % x
            r3test(1,cpt,:) = r3(1,cpt) + sign(nudge)*dnudges;
        case 2
            % y
            r3test(2,cpt,:) = r3(2,cpt) + sign(nudge)*dnudges;
        case 3
            % z
            r3test(3,cpt,:) = r3(3,cpt) + sign(nudge)*dnudges;
        otherwise
            error
    end
    clear r3
    
    im.h.raw = frame.h;
    im.h.s = imfilter(im.h.raw,handles.gaussian);
    [im.h.dx,im.h.dy] = gradient(im.h.s);
    if ~handles.tracking2D
    im.v.raw = frame.v;
    im.v.s = imfilter(im.v.raw,handles.gaussian);
    [im.v.dx,im.v.dy] = gradient(im.v.s);
    else
        im.v=[];
    end
    E = zeros(1,length(dnudges));
    for nidx = 1:length(dnudges)
        r3new = r3test(:,:,nidx);
        z = [0 0 r3new(:,2)' 0 0];
        np0 = normal_plane(r3new,0);   % np0 ~ column vectors that span the plane
        np1 = normal_plane(r3new,1);
        [E(nidx),~] = Bimfun_wm4(z, r3new, np0, np1, im.h, im.v, handles.dt, ...
            handles.trpmtrs.sigma_prior, handles.trpmtrs.sigma2_prior, handles.calib);
    end
    [E,midx] = min(E);
    r3new = r3test(:,:,midx);
    clear dnudges r3test im nidx z np0 np1 midx
    
    handles.whisker(handles.current_whisker).Emin(frameidx) = E;
    titlestring = sprintf('%d:E=%.2f',frameidx,E);
    set(handles.currentframe_display,'string',titlestring)
    clear  titlestring
    r2new = projection2(r3new,handles.calib);
    axes(handles.viewh)
    plot(r2new(1,cpt,1),r2new(2,cpt,1),'.w')
    if ~handles.tracking2D
        axes(handles.viewv)
        plot(r2new(1,cpt,2),r2new(2,cpt,2),'.w')
    end
%     lastframeidx = modsubtract(frameidx,1,handles.video.nframes);
    
    
    % Find intersection between whisker polynomial and folliclemask
    %     roi.x = round(handles.fp_all(handles.frameidx-1,1,cview)) + [-40:40];
    %     roi.y = round(handles.fp_all(handles.frameidx-1,2,cview)) + [-40:40];
    %     roi.x = roi.x((roi.x>=1)&(roi.x<=size(frame,2)));
    %     roi.y = roi.y((roi.y>=1)&(roi.y<=size(frame,1)));
    %     if handles.frameidx == handles.video.startframe
    %         if ~isempty(handles.folliclemask{lastframe,cview})
    %             roi.x = round(handles.whisker(handles.current_whisker).fp_all(1,lastframe,cview)) + [-40:40];
    %             roi.y = round(handles.whisker(handles.current_whisker).fp_all(2,lastframe,cview)) + [-40:40];
    %             %roi.y = round(handles.fp_all(2,lastframe)) + [-40:40];
    %             folliclemask_old = handles.folliclemask{lastframe,cview};
    %             fpidx_old = handles.whisker(handles.current_whisker).fpidx(lastframe,cview);
    %         else
    %             roi.x = [50:size(frame,2)];%50
    %             roi.y = [50:size(frame,1)];%50
    %             folliclemask_old = [];
    %             fpidx_old = [];
    %         end
    %     else
    
    % 210216 following deleted...
    %         roi.x = [50:size(frame,2)];%50
    %         roi.y = [50:size(frame,1)];%50
    %         folliclemask_old = [];
    %         fpidx_old = [];
    % %     end
    %     roi.x = roi.x(roi.x<=size(frame,2));
    %     roi.y = roi.y(roi.y<=size(frame,1));
    %
    %     [p,~] = fit_curve(rnew,handles.polyorder);
    %     theta = base_angle(p,0);
    %     clear p
    %     % 140116: Nb the frameidx-1 in following line fails at startframe
    %     % (should use the modsubtract function)
    %     folliclemask = snout_segment(frame,roi,...
    %         folliclemask_old,fpidx_old,theta,snout_orientation,handles);
    %     clear roi theta
    %
    %     %         [folliclemask{1}] = snout_segment(double(handles.frame.h(:,:,1)),roi.h,folliclemask_old{1}, fpidx_old{1}, theta(1),'horizontal');
    %
    %
    
    % Postprocess:
    %     [~, p, fp, fpidx, sfp, s] = postproc(rnew, folliclemask, handles.trpmtrs.s0, handles.dt, handles.polyorder);
    %     r = rnew;
    %     [rnew, fp, fpidx, sfp, s, sp2] = postproc_v3(rnew, folliclemask, s0, handles.dt);
    
    %     [r3new, fp3, ~, sfp, s, sp2] = ...
    %         postproc_wm4a(r3new, handles.folliclemask{frameidx,1}, handles.trpmtrs.s0, handles.dt);
    
    %             postproc(rnew(:,:,2), folliclemask{2}, handles.trpmtrs.s0.v, handles.dt, handles.trpmtrs.sstar.v);
    
    doplot = true;  % should really be made consistent with 'doplot' in the main loop
    if doplot
        %         b = zeros(2,length(s)+1);
        %         b(1,:) = polyval(p(1,:),[sfp s]);
        %         b(2,:) = polyval(p(2,:),[sfp s]);
        t = [0:handles.dt:1];
        b3 = bezierval(r3new,t);
%         b2 = zeros(2,size(b3,2),2);
%         b2(:,:,1) = b3(1:2,:);
%         b2(:,:,2) = handles.calib.matrix*b3 + handles.calib.vector*ones(1,size(b3,2));
%         r2 = zeros(2,3,2);
%         r2(:,:,1) = r3new(1:2,:);
%         r2(:,:,2) = handles.calib.matrix*r3new + handles.calib.vector*ones(1,3);
        r2 = projection2(r3new,handles.calib);
        b2 = projection2(b3,handles.calib);
        clear t
        axes(handles.viewh)
        plot(b2(1,:,1),b2(2,:,1),'r',r2(1,:,1),r2(2,:,1),'r.')%,fp2(1,1),fp2(2,1),'y.')
        plot(r2new(1,:,1),r2new(2,:,1),'msq'),
        if ~handles.tracking2D
        axes(handles.viewv)
        plot(b2(1,:,2),b2(2,:,2),'r',r2(1,:,2),r2(2,:,2),'r.')%,fp2(1,2),fp2(2,2),'y.')
        plot(r2new(1,:,2),r2new(2,:,2),'msq'),
        end
        %         plot(polyval(p(1,:),handles.trpmtrs.sstar),polyval(p(2,:),handles.trpmtrs.sstar),'go')
        clear b3
    end
    clear sfp s doplot
    
    handles.whisker(handles.current_whisker).r3all(frameidx,:,:) = r3new;
    %     handles.whisker(handles.current_whisker).fp_all(handles.frameidx,:,cview) = fp;
    %     handles.folliclemask{handles.frameidx,cview} = folliclemask;
    %     handles.whisker(handles.current_whisker).fpidx_all(:,handles.frameidx,cview) = fpidx;
    handles.whisker(handles.current_whisker).tracked(frameidx) = true;
    %     tstar = 0;
    %     handles.whisker(handles.current_whisker).kappa_all(handles.frameidx,cview) = curvature(rnew,tstar);
    %     handles.whisker(handles.current_whisker).theta_all(handles.frameidx,cview) = base_angle(rnew,tstar);
    clear fp fpidx
    
    % Save to disk:
    whisker = handles.whisker;
    %     rall = handles.rall;
    %     fp_all = handles.fp_all;
    %     kappa_all = handles.kappa_all;
    %     theta_all = handles.theta_all;
    trpmtrs = handles.trpmtrs;
    roi = handles.roi;
    calib = handles.calib;
    %     tracked = handles.tracked;
    %     pmtrs.sstar = handles.trpmtrs.s0;
    %     pmtrs.sigma_prior = handles.trpmtrs.sigma_prior;
    %     pmtrs.subtract_image_mean = handles.trpmtrs.subtract_image_mean;
    %     save(handles.trfname,'rall','kappa_all','theta_all','fp_all','trpmtrs','tracked')
    save(handles.trfname,'whisker','trpmtrs','roi','calib')
    clear rall kappa_all theta_all fp_all trpmtrs roi calib
    
end


guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_move_control_point.
function pushbutton_move_control_point_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_move_control_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

error('Update the curve to postprocess contour')

handles = move_control_point_Callback(handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = move_control_point_Callback(handles)

cpt = handles.current_pt;
if cpt==0
    return
end

axes(handles.viewv)
cla
imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1)),
colormap gray
title(handles.frameidx)
hold on

r = squeeze(handles.rall(handles.frameidx,:,:));
plot(r(1,:),r(2,:),'r.-')
plot(r(1,cpt),r(2,cpt),'wo')
[x,y] = ginput(1);
plot(x,y,'w+')
r(:,cpt) = [x y]';
plot(r(1,:),r(2,:),'r--o')

D = size(r,2);
handles.rall(handles.frameidx,:,:) = r;

clear cpt r ddr

% guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in radiobutton_kinematics_test.
function radiobutton_kinematics_test_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_kinematics_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_kinematics_test

% --- Executes on button press in radiobutton_continuous_tracking.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton_continuous_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_continuous_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_continuous_tracking

% dummy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in radiobutton_subtract_image_mean.
function radiobutton_subtract_image_mean_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_subtract_image_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_subtract_image_mean

handles.meanframe = init_subtract_image_mean(get(hObject,'value'),handles.roi,handles.video);

% meanframe.h = zeros(handles.video.height.h,handles.video.width.h,3);
% meanframe.v = zeros(handles.video.height.v,handles.video.width.v,3);
%
% switch get(hObject,'Value')
%     case true
%         % subtract mean image
%         set(handles.text_status,'BackgroundColor',handles.red)
%         set(handles.text_status,'String','BUSY')
%         pause(.01)
%         space = 10;
%         for frameidx = 1:space:handles.video.nframes
%             tmp = load_frame_2views(handles.video, frameidx, handles.roi);
%             meanframe.h = meanframe.h + double(tmp.h)/(handles.video.nframes/space);
%             meanframe.v = meanframe.v + double(tmp.v)/(handles.video.nframes/space);
%         end
%         clear space frameidx
%         set(handles.text_status,'BackgroundColor',handles.green)
%         set(handles.text_status,'String','idle')
%     case false
%         % do not subtract mean image
%     otherwise
%         error('This should never happen')
% end
%
% fn = [handles.fname(1:end-3) 'meanframe'];
% save(fn,'meanframe')
% handles.meanframe = meanframe;
% clear meanframe fn

titles.h = ''; titles.v = '';
% plot_currentframe(handles,titles)
plot_frames(frame,handles.video,handles.meanframe,titles,handles)


% Update the energy threshold:
for w = 1:handles.Nwhiskers
    if get(handles.radiobutton_subtract_image_mean,'Value')
        handles.whisker(w).energy_threshold = handles.default_energy_threshold_subtract_image_mean;
    else
        handles.whisker(w).energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
    end
end
set(handles.edit_current_whisker_energy_threshold,'String', handles.whisker(handles.current_whisker).energy_threshold);

guidata(hObject, handles);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanframe = init_subtract_image_mean(radiobutton_value,roi,video)

width.h = roi.h(2)-roi.h(1)+1;
height.h = roi.h(4)-roi.h(3)+1;
meanframe.h = zeros(height.h,width.h,3);

if isfield(roi, 'v')
width.v = roi.v(2)-roi.v(1)+1;
height.v = roi.v(4)-roi.v(3)+1;
meanframe.v = zeros(height.v,width.v,3);
end
clear height width

switch radiobutton_value
    case true
        % subtract mean image
        set(handles.text_status,'BackgroundColor',handles.red)
        set(handles.text_status,'String','BUSY')
        pause(.01)
        space = 10;
        %         frameidx = 1;   % initialisation
        %         tmp = load_frame_2views(handles.video, frameidx, handles.roi,get(handles.checkbox_rotatevview,'value'));
        %         meanframe.h = double(tmp.h)/(handles.video.nframes/space);
        %         meanframe.v = double(tmp.v)/(handles.video.nframes/space);
        for frameidx = 1:space:video.h.nframes
            tmp = load_frame(video.h, frameidx, roi.h);%,get(handles.checkbox_rotatevview,'value'));
            meanframe.h = meanframe.h + double(tmp)/(video.h.nframes/space);
        end
        clear frameidx tmp
        for frameidx = 1:space:video.v.nframes
            tmp = load_frame(video.v, frameidx, roi.v);%,get(handles.checkbox_rotatevview,'value'));
            meanframe.v = meanframe.v + double(tmp)/(video.v.nframes/space);
        end
        clear frameidx tmp
        clear space frameidx
        set(handles.text_status,'BackgroundColor',handles.green)
        set(handles.text_status,'String','idle')
    case false
        % do not subtract mean image
    otherwise
        error('This should never happen')
end

% handles.meanframe = meanframe;
% clear meanframe fn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_frames(frame,video,meanframe,titles,haxes)

axes(haxes.h)
cla
% if video. v exist
if isfield(video,'v')
    xlim([1 max([video.h.width.roi,video.v.width.roi])])
    ylim([1 max([video.h.height.roi,video.v.height.roi])])
else
    xlim([1 video.h.width.roi])
    ylim([1 video.h.height.roi])
end
set(gca,'ydir','reverse')
axis image off
hold on
imagesc(double(frame.h(:,:,1))-meanframe.h(:,:,1));
title(titles.h)
colormap gray

if isfield(video,'v')
axes(haxes.v)
cla
xlim([1 max([video.h.width.roi,video.v.width.roi])])
ylim([1 max([video.h.height.roi,video.v.height.roi])])
set(gca,'ydir','reverse')
axis image off
hold on
imagesc(double(frame.v(:,:,1))-meanframe.v(:,:,1));
title(titles.v)
colormap gray
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_image(image,handles)

axes(handles.viewh)
cla
xlim([1 max([handles.hvideo.width.raw,handles.vvideo.width.raw])])
ylim([1 max([handles.hvideo.height.raw,handles.vvideo.height.raw])])
set(gca,'ydir','reverse')
axis image off
hold on
imagesc(image.h.s);
colormap gray

axes(handles.viewv)
cla
xlim([1 max([handles.hvideo.width.raw,handles.vvideo.width.raw])])
ylim([1 max([handles.hvideo.height.raw,handles.vvideo.height.raw])])
set(gca,'ydir','reverse')
axis image off
hold on

% if get(handles.checkbox_rotatevview,'value')
%     imagesc(image.v.s');
% else
imagesc(image.v.s);
% end
colormap gray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.Nwhiskers = 1;
% handles.whisker(1,1:handles.Nwhiskers).label = 'xx';
% if get(handles.radiobutton_subtract_image_mean,'Value')
%     handles.whisker(1,1:handles.Nwhiskers).energy_threshold = handles.default_energy_threshold_subtract_image_mean;
% else
%     handles.whisker(1,1:handles.Nwhiskers).energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
% end
% handles.whisker(1,1:handles.Nwhiskers).selected = 1;
% handles.whisker(1,1:handles.Nwhiskers).tracked = zeros(1,handles.video.nframes);
% handles.whisker(1,1:handles.Nwhiskers).r3all = zeros(handles.video.nframes,3,handles.n_control_pts);
% handles.whisker(1,1:handles.Nwhiskers).fp3_all = zeros(handles.video.nframes,3);
% handles.whisker(1,1:handles.Nwhiskers).s0 = [];

% whisker = handles.whisker;
% trpmtrs = handles.trpmtrs;
% roi = handles.roi;
% calib = handles.calib;
% save(handles.trfname,'-mat','whisker','trpmtrs','roi','calib')
% clear trpmtrs whisker roi calib

handles.frameidx = handles.video.startframe;  % current frame
handles.frame = load_frame_2views(handles.video,handles.frameidx,handles.roi);%,get(handles.checkbox_rotatevview,'value'));
set(handles.currentframe_display,'string',num2str(handles.frameidx))
handles.meanframe = init_subtract_image_mean(get(handles.radiobutton_subtract_image_mean,'value'),handles.roi,handles.video);
titles.h = ''; titles.v = '';
% plot_currentframe(handles, titles);
plot_frames(handles.frame,handles.video,handles.meanframe,titles,handles)

delete(handles.trfname)
handles = initialise_new_track(handles);

% handles.trfname
%
% handles.rall = zeros(handles.video.nframes,2,handles.n_control_pts);
% handles.kappa_all = zeros(1,handles.video.nframes);
% handles.theta_all = zeros(1,handles.video.nframes);
% handles.fp_all = zeros(2,handles.video.nframes);
% handles.fpidx_all = zeros(1,handles.video.nframes);
% handles.Emin = zeros(1,handles.video.nframes);
% handles.tracked = zeros(1,handles.video.nframes);
%
% handles.current_pt = 0;
% handles.trpmtrs.s0 = [];
%
% handles.frameidx = handles.video.startframe;
% handles.frame = load_frame(handles.video,handles.frameidx);
% set(handles.currentframe_display,'string',num2str(handles.frameidx))
%
% axes(handles.viewh)
% cla
% image(handles.frame)
%
% axes(handles.viewv)
% cla
% imagesc(double(handles.frame(:,:,1))-handles.meanframe(:,:,1))
% % title(frameidx)
% colormap gray
% hold on
%
handles.continuous_tracking = false;
set(handles.radiobutton_continuous_tracking,'value',handles.continuous_tracking)
% handles.trpmtrs.subtract_image_mean = true;
% set(handles.radiobutton_subtract_image_mean,'value',handles.trpmtrs.subtract_image_mean)

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_plot_kinematics.
function pushbutton_plot_kinematics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_kinematics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% rall = handles.rall;
%
% nframes = size(rall,1);
frames = 1:handles.video.h.nframes;

% D = size(rall,3);

% idx
% size(handles.tracked)
% size(handles.theta_all)

% % compute base angle
% theta = zeros(length(idx),D-1);
% for d = 1:D-1
%     dy = rall(idx,2,d+1)-rall(idx,2,1);
%     dx = rall(idx,1,d+1)-rall(idx,1,1);
%     theta(idx,d) = atan2(-dy,dx)*(180/pi);
% end
% clear d dx dy

% % compute curvature: kappa=(x^'y^('')-y^'x^(''))/((x^('2)+y^('2))^(3/2)).
% tmp = diff(rall,1,3);
% dt = reshape(sum(tmp.^2,2),nframes,D-1);    % differential lengths along the contours
% clear tmp
% % t = cumsum(dt,2);   % distance along the contours
% dxdt = reshape(diff(rall(:,1,:),1,3),nframes,D-1) ./ dt;
% dxdt = dxdt(:,2:end);
% d2xdt2 = reshape(diff(rall(:,1,:),2,3),nframes,D-2) ./ dt(:,2:end).^2;
% dydt = reshape(diff(rall(:,2,:),1,3),nframes,D-1) ./ dt;
% dydt = dydt(:,2:end);
% d2ydt2 = reshape(diff(rall(:,2,:),2,3),nframes,D-2) ./ dt(:,2:end).^2;
% kappa = dxdt.*d2ydt2 - dydt.*d2xdt2;
% kappa = kappa ./ (dxdt.^2+dydt.^2).^(3/2);

f1=figure;
% if ~handles.tracking2D
 spr = 2;
% else
% spr=2;
% end
handles.legend={};
for w = 1:handles.Nwhiskers
    if ~handles.whisker(w).selected
        continue
    end
    handles.legend=[handles.legend {handles.whisker(w).label}];
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    idx = find(handles.whisker(w).tracked);
    %     theta1 = handles.whisker(w).theta_all(idx,1);
    %     theta2 = handles.whisker(w).theta_all(idx,2);
    %     theta1(theta1<0) = theta1(theta1<0)+360;
    %     theta2(theta2<0) = theta2(theta2<0)+360;
    %     kappa1 = handles.whisker(w).kappa_all(idx,1);
    %     kappa2 = handles.whisker(w).kappa_all(idx,2);
    theta = zeros(3,length(handles.whisker(w).tracked));
    azimuth = zeros(1,length(handles.whisker(w).tracked));
    elevation = zeros(1,length(handles.whisker(w).tracked));
    twist = zeros(1,length(handles.whisker(w).tracked));
    curv_hc = zeros(3,length(handles.whisker(w).tracked));  % head-centred
    curv_fc = zeros(3,length(handles.whisker(w).tracked));  % follicle-centred - later
    curv3 = zeros(1,length(handles.whisker(w).tracked));
    np = zeros(2,length(handles.whisker(w).tracked));
    r3 = squeeze(handles.whisker(w).r3all);
    for i = 1:length(idx)
        fr = idx(i);
        r = squeeze(handles.whisker(w).r3all(fr,:,:));
        theta(:,fr) = base_angle3(r,0);
        curv_hc(1,fr) = curvature(r([2 3],:),0);
        curv_hc(2,fr) = curvature(r([1 3],:),0);
        curv_hc(3,fr) = curvature(r([1 2],:),0);
        curv3(fr) = curvature3(r,0);
        tv = bezierdtval(r,0);  % tangent to whisker at base
        tv = tv/norm(tv);
        cv = bezierdt2val(r,0); % 2nd derivative to whisker at base
        NPO = eye(3)-tv*tv'; % projection operator, onto plane normal to 'tv'
        %         xpr = NPO*[1 0 0]'; % projection of x axis
        %         zpr = NPO*[0 0 1]'; % projection of z axis
        np(1,fr) = [1 0 0]*NPO*cv;   % projection of cv onto NP, projected onto x axis
        np(2,fr) = [0 0 1]*NPO*cv;   % projection of cv onto NP, projected onto z axis
        %         curv(1,fr) = cv'*xpr; % projection of 2nd deriv vector
        %         curv(2,fr) = cv'*zpr;
    end
    clear i fr tv cv
    
    azimuth(idx) = theta(3,idx);
    elevation(idx) = theta(1,idx);
    twist(idx) = atan2(-np(2,idx),np(1,idx))*(180/pi);
    if handles.tracking2D
        figure(f1)
        h1(1) = subplot(spr,1,1); hold on
        plot(frames(idx),azimuth(idx),handles.colours{wmod}), ylabel('Azimuth angle [\circ]'),legend(handles.legend)
        title('Azimuth: 90\circ is normal to ant-post axis; >90 \circ for whisker protracted')
        h1(2) = subplot(spr,1,2); hold on
        plot(frames(idx),curv_hc(3,idx),handles.colours{wmod}), ylabel('\kappa Horizontal [1/pixel]'),xlabel('Frame number')
        linkaxes(h1,'x')


    else
        figure(f1)
        h1(1) = subplot(spr,2,1); hold on
        plot(frames(idx),azimuth(idx),handles.colours{wmod}), ylabel('Azimuth angle [\circ]'),legend(handles.legend)
        title('Azimuth: 90\circ is normal to ant-post axis; >90\circ for whisker protracted')
        h1(2) = subplot(spr,2,3); hold on
        plot(frames(idx),elevation(idx),handles.colours{wmod}), ylabel('Elevation angle [\circ]'),xlabel('Frame number')
        title('Elevation: 90\circ is horizontal; >90\circ for whisker tip oriented up')
%         h1(3) = subplot(spr,2,5); hold on
%         plot(frames(idx),twist(idx),handles.colours{wmod}), ylabel('twist')
%         title('twist: 90deg is concave down; <90deg for concave posterior')
        
        h2(1) = subplot(spr,2,4); hold on
        plot(frames(idx),curv_hc(1,idx),handles.colours{wmod}), ylabel('\kappa Coronal [1/pixel]')
        h2(2) = subplot(spr,2,2); hold on
        plot(frames(idx),curv_hc(3,idx),handles.colours{wmod}), ylabel('\kappa Horizontal [1/pixel]'),xlabel('Frame number')
%         h2(3) = subplot(spr,2,6); hold on
%         plot(frames(idx),curv3(idx),handles.colours{wmod}), ylabel('kappa3')
        
        %     h(3) = subplot(spr,1,3); hold on
        %     % plot(frames(idx),medfilt1(kappa(idx,:),5))
        %     % plot(frames(idx),kappa(idx,:),'k',frames(idx),handles.kappa_all(idx),'r')
        %     plot(frames(idx),kappa1,colours{w}), ylabel('h: kappa')
        %     h(4) = subplot(spr,1,4); hold on
        %     plot(frames(idx),kappa2,colours{w}), ylabel('v: kappa')
        % h(3) = subplot(spr,1,3);
        % plot(frames(idx),handles.fp_all(1,idx),'.')
        % title('Follicle X position'), xlabel('Frame'), ylabel('pixels')
        % h(4) = subplot(spr,1,4);
        % plot(frames(idx),handles.fp_all(2,idx),'.')
        % title('Follicle Y position'), xlabel('Frame'), ylabel('pixels')
        %     h(3) = subplot(spr,1,3);
        %     feather(curv(1,idx),curv(2,idx),handles.colours{w})
        %     linkaxes(h,'x')
        
        a = axis;
        tmp = curv_hc(:,idx);
        tmp = [tmp(:)' curv3(idx)];
        %     axis([a(1:2) min(tmp) max(tmp)])
        clear tmp
        linkaxes(h2,'x')
        linkaxes(h1,'x')
        
        
    end
end
clear w


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_sigma_prior_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma_prior as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma_prior as a double

handles.trpmtrs.sigma_prior = str2double(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_sigma_prior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_frame_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frame_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_frame_interval as a double

handles.frame_interval = str2double(get(hObject,'String'));
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function edit_frame_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_auto_initialise.
function radiobutton_auto_initialise_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_auto_initialise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_auto_initialise


% --- Executes on selection change in popupmenu_choose_auto_initialise_file.
function popupmenu_choose_auto_initialise_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_auto_initialise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_choose_auto_initialise_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_choose_auto_initialise_file


% --- Executes during object creation, after setting all properties.
function popupmenu_choose_auto_initialise_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_auto_initialise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.tr4');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose init file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);

clear ff nr_files string_list i



% --- Executes on button press in pushbutton_start_batch.
function pushbutton_start_batch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.radiobutton_auto_initialise,'Value',true)
set(handles.radiobutton_continuous_tracking,'Value',true);
menu = cellstr(get(handles.popupmenu_choose_batch_file,'String'));
batchfn = menu{get(handles.popupmenu_choose_batch_file,'Value')};
batch_fid = fopen(batchfn,'rt');
clear menu
reportfn = 'report.rep';
report_fid = fopen(reportfn,'at');
time = fix(clock);
time = sprintf('%d:%d:%d',time(4:6));
fprintf(report_fid,'\nStarting batch %s at %s on %s...\n',batchfn,time,date);
clear time batch_fn
stopbatch = false;
batch_counter=1;
while (batch_fid && ~stopbatch)
    nextvideo=fgetl(batch_fid);
 
    if nextvideo==-1
        stopbatch = true;
        disp('Batch finished')
        continue
    else
        disp(['Starting video number ' num2str(batch_counter)])
        datfiles=strsplit(nextvideo);
        handles.fname.h =datfiles{1};
        if ~handles.tracking2D
            handles.fname.v =datfiles{2};
            handles.trfname = [handles.fname.h(1:end-4) '.tr4'];
        else
            handles.trfname = [handles.fname.h(1:end-4) '.tr4_2D'];
        end
        %fprintf(report_fid,' Opening video %s\n',handles.fname.h)    % ditto
        
        clear datenum idx ff fn
        clear ff nr_files string_list i
        %    else
        % there are no suitable tr files in current directory.  use default.
        %    end
        %%Initialise horizontal video
        handles.video.h = initialise_new_video(handles.fname.h,handles,handles.dir.h);
        % Set default ROIs:
        handles.roi.h = [1,handles.video.h.width.raw,1,handles.video.h.height.raw];
        handles.video.h.width.roi = handles.roi.h(2)-handles.roi.h(1)+1;
        handles.video.h.height.roi = handles.roi.h(4)-handles.roi.h(3)+1;
        handles.currentframe.h = handles.video.h.startframe;
        handles.mastervideo_selected = true;
        
        if ~handles.tracking2D
            %%Initialise vertical video
            % get a FID and header information ('video' structure):
            handles.video.v = initialise_new_video(handles.fname.v,handles,handles.dir.v);
            handles.roi.v = [1,handles.video.v.width.raw,1,handles.video.v.height.raw];
            handles.video.v.width.roi = handles.roi.v(2)-handles.roi.v(1)+1;
            handles.video.v.height.roi = handles.roi.v(4)-handles.roi.v(3)+1;
            handles.currentframe.v = handles.video.v.startframe;
        end
        % user should already have selected the "master" horizontal view video
        if handles.mastervideo_selected
            if handles.tracking2D
                % make the master frameidx table
                nframes = handles.video.h.nframes;
                frameidxs.h = modadd(1:nframes,handles.video.h.startframe-1,nframes);
                idx = find(frameidxs.h==1);
                clear frameidxs idx
            else
                % make the master frameidx table
                if handles.video.h.nframes==handles.video.v.nframes
                    nframes = handles.video.h.nframes;
                    % straightforward
                    frameidxs.h = modadd(1:nframes,handles.video.h.startframe-1,nframes);
                    frameidxs.v = modadd(1:nframes,handles.video.v.startframe-1,nframes);
                    idx = find(frameidxs.h==1);
                    handles.framenums_vfromh = modadd(1:nframes,frameidxs.v(idx)-1,nframes);
                    clear frameidxs idx
                else
                    error('Extend the code!')
                end
            end
            % initialise parameters that apply at the level of the 2 views together
            handles.dt = 0.01;
            hsize = [3 3]; % default
            sigma = .05;
            handles.gaussian = fspecial('gaussian', hsize, sigma);
            clear hsize sigma
            % meanframe
            meanframe = init_subtract_image_mean(false,handles.roi,handles.video);
            % 060416: is following line needed?
            fn = [handles.fname.h(1:end-3) 'meanframe'];    % NB H is the master
            save(fn,'meanframe')
            handles.meanframe = meanframe;
            clear fn meanframe
            % load and plot the first frame:
            titles.h = sprintf('Frame %d',handles.currentframe.h);
            haxes.h = handles.viewh;
            if ~handles.tracking2D
                titles.v = sprintf('Frame %d',handles.currentframe.v);
                haxes.v = handles.viewv;
            end
            handles.frame = load_and_plot_frames(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,haxes);
            clear haxes
            % set default tracking file
            %handles.trfname = [handles.fname.h(1:end-4) '.tr4'];
            % If a 'tr' file for current video exists, load the data; else, initialise:
            handles = initialise_new_track(handles);
            handles.current_view = 1;
            handles.current_pt = 0;
            if isfield(handles,'calib')
                plot_contours(handles.currentframe,handles)
            end
        else
            error('Select the horizontal view (master) video first')
        end
        %%fit snakes
        set(handles.pushbutton_fit_snakes,'enable','on')
        [handles,fnoutput] = fit_snakes(handles);
        switch fnoutput
            case -1
                % should only happen if user switched off the
                % continuous_tracking radio button
                stopbatch = true;
                fprintf(report_fid,' Track failure: output=%d\n',fnoutput);
            case 1
                % fit_snakes successfully reached end of video
                disp('Tracker reached end of video without error')
                fprintf(report_fid,' Tracker reached end of video without error\n');
            case 2
                % Tracking failed somewhere in the video
                %stopbatch = true;
                fprintf(report_fid,' Track failure: output=%d\n',fnoutput);
            otherwise
                error('Unrecognised output from fit_snakes %d', fnoutput)
        end
        if ~get(handles.radiobutton_continuous_tracking,'Value')
            % ie user interrupted the track
            stopbatch = true;
            fprintf(report_fid,' Track failure: user abort\n');
        end
    end
    batch_counter=batch_counter+1;
end
fclose(batch_fid);
fclose(report_fid);
clear batch_fid report_fid
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_choose_batch_file.
function popupmenu_choose_batch_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_batch_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_choose_batch_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_choose_batch_file

ff = dir('*.bat');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose batch file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i


% --- Executes during object creation, after setting all properties.
function popupmenu_choose_batch_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_choose_batch_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.bat');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose batch file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function edit_miumax_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_miumax (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_miumax as text
% %        str2double(get(hObject,'String')) returns contents of edit_miumax as a double
% 
% handles.trpmtrs.miumax = str2double(get(hObject,'String'));
% guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
% function edit_miumax_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_miumax (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function z = modsubtract(x,y,n)
% modular subtraction.
% eg if x=5,y=2,n=9, z = 3 (as normal)
% eg  if x=2,y=2,n=9, z = 9 (wrap-around)

z = x-y;
% if z <= 0
%     z = z+n;
% end
idx =  z <= 0;
z(idx) = z(idx) + n;
clear idx


function z = modadd(x,y,n)
% modular addition
% eg if x=5,y=2,n=9, z = 7 (as normal)
% eg  if x=8,y=2,n=9, z = 1 (wrap-around)

z = x+y;

idx =  z > n;
z(idx) = z(idx) - n;
clear idx



function handles = initialise_new_track(handles)
% if tracking data file exists, load both data and pmtrs:
if exist(handles.trfname,'file')
    %     load(handles.trfname,'rall','kappa_all','theta_all','fp_all','trpmtrs','tracked','-mat')
    load(handles.trfname,'whisker','trpmtrs','roi','calib','-mat')
    handles.whisker = whisker;
    handles.roi = roi;
    if exist('calib','var')
        handles.calib = calib;
    end
    clear whisker roi calib

    handles.n_control_pts = size(handles.whisker(1).r3all,3);
    if handles.n_control_pts~=3
        error('Bad tr file: N control points must be 3')
    end
    handles.Nwhiskers = length(handles.whisker);
    
    if ~isfield(trpmtrs,'tracking_direction')   % new in v4b2 050916
        % for backwards compatibility with tr files that might lack this pmtr
        trpmtrs.tracking_direction = handles.trpmtrs.tracking_direction;
    end
    
    handles.trpmtrs = trpmtrs;
    set(handles.edit_sigma_prior,'String',handles.trpmtrs.sigma_prior);
    if isfield(handles,'trpmtrs_sigma2_prior')
        set(handles.edit_sigma2_prior,'String',handles.trpmtrs_sigma2_prior);
    end
    if isfield(handles,'trpmtrs_snout_sigma')
        set(handles.edit_snout_sigma,'String',handles.trpmtrs_snout_sigma);
    end
    %set(handles.radiobutton_subtract_image_mean,'value',handles.trpmtrs.subtract_image_mean)
    
    set(handles.pushbutton_tracking_direction,'String',handles.trpmtrs.tracking_direction)  % new in 4b2 (050916)

%     set(handles.edit_energy_threshold,'String',handles.trpmtrs.energy_threshold);
%     if isfield(trpmtrs,'miumax')
%         set(handles.edit_miumax,'String',handles.trpmtrs.miumax);
%     else
%         handles.trpmtrs.miumax = str2double(get(handles.edit_miumax,'String'));
%     end
    clear trpmtrs
else    % if no tracking data file exists, initialise.
    handles.n_control_pts = 3;
    handles.Nwhiskers = 1;

    % use default roi
%     roi.h = handles.roi.h;
%     roi.v = handles.roi.v;
    roi = handles.roi;
%     handles.roi = roi;

    % whisker-level pmtrs
    whisker(1,1:handles.Nwhiskers).label = 'xx';
    whisker(1,1:handles.Nwhiskers).energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;

    whisker(1,1:handles.Nwhiskers).selected = 1;
    whisker(1,1:handles.Nwhiskers).tracked =  zeros(1,handles.video.h.nframes);
    whisker(1,1:handles.Nwhiskers).r3all = zeros(handles.video.h.nframes,3,handles.n_control_pts);
    %     whisker(1,1:handles.Nwhiskers).kappa_all = zeros(handles.video.nframes,2);
    %     whisker(1,1:handles.Nwhiskers).theta_all = zeros(handles.video.nframes,2);
    whisker(1,1:handles.Nwhiskers).fp3_all = zeros(handles.video.h.nframes,3);
    whisker(1,1:handles.Nwhiskers).s0 = [];
    handles.whisker = whisker;

    % video-level pmtrs
    trpmtrs.sigma_prior = str2double(get(handles.edit_sigma_prior,'String'));
    trpmtrs.sigma2_prior = str2double(get(handles.edit_sigma2_prior,'String'));
    trpmtrs.snout_sigma = str2double(get(handles.edit_snout_sigma,'String'));
    trpmtrs.tracking_direction = get(handles.pushbutton_tracking_direction,'String'); % new in 4b2 (050916)

%     trpmtrs.energy_threshold = str2double(get(handles.edit_energy_threshold,'String'));
%     trpmtrs.miumax = str2double(get(handles.edit_miumax,'String'));
    trpmtrs.subtract_image_mean = false;
    handles.trpmtrs = trpmtrs;
    
    save(handles.trfname,'-mat','whisker','trpmtrs','roi')
    clear trpmtrs whisker roi
end

% handles.video.width.h = handles.roi.h(2)-handles.roi.h(1)+1;
% handles.video.width.v = handles.roi.v(2)-handles.roi.v(1)+1;
% handles.video.height.h = handles.roi.h(4)-handles.roi.h(3)+1;
% handles.video.height.v = handles.roi.v(4)-handles.roi.v(3)+1;

if isfield(handles,'calib')
    set(handles.pushbutton_fit_snakes,'enable','on')
end

set(handles.text_Nwhiskers,'String',handles.Nwhiskers)
handles.current_whisker = 1;
set(handles.text_current_whisker,'String',handles.current_whisker)
set(handles.edit_current_whisker_label,'String',handles.whisker(handles.current_whisker).label)
set(handles.radiobutton_current_whisker_selected,'value',handles.whisker(handles.current_whisker).selected)
set(handles.edit_current_whisker_energy_threshold,'String', handles.whisker(handles.current_whisker).energy_threshold);




% 180216: not sure following is used; try commenting out:
% handles.polyorder = handles.n_control_pts-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frame] = load_and_plot_frames(video,framenum,roi,meanframe,titles,haxes)

frame.h = load_frame(video.h,framenum.h,roi.h);%,get(handles.checkbox_rotatevview,'value'));

if isfield(video,'v')
    frame.v = load_frame(video.v,framenum.v,roi.v);
end
% set(handles.currentframe_display,'string',num2str(frameidx))

% video.width = size(frame,2);
% % video.width.v = size(frame.v,2);
% video.height = size(frame,1);
% video.height.v = size(frame.v,1);

% apply currently selected 'subtract image mean' setting:
% meanframe = init_subtract_image_mean(get(handles.radiobutton_subtract_image_mean,'value'),handles);

plot_frames(frame,video,meanframe,titles,haxes)



function edit_specify_tr_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_specify_tr_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_specify_tr_file as text
%        str2double(get(hObject,'String')) returns contents of edit_specify_tr_file as a double

handles.trfname = get(hObject,'String');
handles = initialise_new_track(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_specify_tr_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_specify_tr_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_select_rois.
function pushbutton_select_rois_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUID ATA)

debug = 0;

roi_old = handles.roi;

if debug
    tmp.h = [300 600 1 400];
    tmp.v = [1 300 1 400];
else
    tmp = setroi_gui(handles.frame.raw,handles.roi);
end
% NB these are coordinates referenced to the raw image

handles.roi = tmp;
roi = tmp;
save(handles.trfname,'roi','-append')
clear tmp roi

% pmtrs that depend on roi need to be reset:
% handles.video.width.h = handles.roi.h(2)-handles.roi.h(1)+1;
% handles.video.width.v = handles.roi.v(2)-handles.roi.v(1)+1;
% handles.video.height.h = handles.roi.h(4)-handles.roi.h(3)+1;
% handles.video.height.v = handles.roi.v(4)-handles.roi.v(3)+1;

handles.meanframe = init_subtract_image_mean(get(handles.radiobutton_subtract_image_mean,'value'),handles.roi,handles.video);

% load and plot the first frame:
titles.h = ''; titles.v = '';
handles.frame = load_and_plot_frame(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,handles);
clear titles

% translate the coordinate system of the tracking solution to that of the
% new roi:
handles.whisker = translate_contours_to_new_roi(handles.whisker,handles.roi, roi_old);

% Since ROIs may have changed, re-compute the 3D-to-vview projection:
if isfield(handles, 'calib')
    handles.calib = fit_projection_from_calibration_data(handles.roi,handles.calib.C);
end

% Replot the tracking solution (if any)
if isfield(handles, 'calib')
    plot_contours(handles.currentframe,handles)
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calib = fit_projection_from_calibration_data(roi,C)

% C is in raw image coords.  Translate xy and XY components to coord frame of
% the xy ROI, and vw components to coord frame of the vw ROI:
origin_XY = roi.h([1,3]);
origin_xy = roi.h([1,3]);
origin_vw = roi.v([1,3]);
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

% Do the regression for 'w':
b = regress(w,[P ones(n,1)]);
calib.mw = b(1:3)';
calib.ow = b(4);
clear b

calib.matrix = [calib.mv;calib.mw];
calib.vector = [calib.ov;calib.ow];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function whisker_new = translate_contours_to_new_roi(whisker_old,roi_new,roi_old)

whisker_new = whisker_old;
Nwhiskers = length(whisker_old);

for w = 1:Nwhiskers
    r3all = whisker_old(w).r3all;
    fp3_all = whisker_old(w).fp3_all;
    tridx = find(whisker_old(w).tracked);
    r3all(tridx,1,:) = r3all(tridx,1,:) + roi_old.h(1) - roi_new.h(1);
    r3all(tridx,2,:) = r3all(tridx,2,:) + roi_old.h(3) - roi_new.h(3);
    fp3_all(tridx,1) = fp3_all(tridx,1) + roi_old.h(1) - roi_new.h(1);
    fp3_all(tridx,2) = fp3_all(tridx,2) + roi_old.h(3) - roi_new.h(3);
    
    %     % For v view, take coordinate rotation into account:
    %     rall.v(tridx,1,:) = rall.v(tridx,1,:) + roi_old.v(3) - roi_new.v(3);
    %     rall.v(tridx,2,:) = rall.v(tridx,2,:) + roi_old.v(1) - roi_new.v(1);
    %     fp_all.v(tridx,1) = fp_all.v(tridx,1) + roi_old.v(3) - roi_new.v(3);
    %     fp_all.v(tridx,2) = fp_all.v(tridx,2) + roi_old.v(1) - roi_new.v(1);
    ! z coord should also be adjusted, but leave this for now
    
    whisker_new(w).r3all(:,:,:,1) = r3all;
    whisker_new(w).fp3_all(:,:,1) = fp3_all;
    clear r3all fp3_all tridx
end


% --- Executes on mouse press over axes background.
function viewh_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to viewh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_rotatevview.
% function checkbox_rotatevview_Callback(hObject, eventdata, handles)
% % hObject    handle to checkbox_rotatevview (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hint: get(hObject,'Value') returns toggle state of checkbox_rotatevview
%
% % If user toggles this halfway through a track, it will create chaos, so
% % don't let them:
% set(hObject,'Value','true')

% 
% function edit_energy_threshold_v_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_energy_threshold_v (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_energy_threshold_v as text
% %        str2double(get(hObject,'String')) returns contents of edit_energy_threshold_v as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit_energy_threshold_v_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_energy_threshold_v (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



% function edit_miumax_v_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_miumax_v (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_miumax_v as text
% %        str2double(get(hObject,'String')) returns contents of edit_miumax_v as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit_miumax_v_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_miumax_v (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function edit_snout_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_snout_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_snout_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_snout_sigma as a double


% --- Executes during object creation, after setting all properties.
function edit_snout_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_snout_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sigma2_prior_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma2_prior_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma2_prior as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma2_prior as a double

handles.trpmtrs.sigma2_prior = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_sigma2_prior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma2_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_play_video.
function pushbutton_play_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   ecs

for w = 1:handles.Nwhiskers
    if ~handles.whisker(w).selected
        break
    end
    clf
    tridx{w} = find(handles.whisker(w).tracked);
    r3{w} = handles.whisker(w).r3all;
    tmp = r3{w}(tridx{w},1,:);
    xmin(w) = min(tmp(:)); xmax(w) = max(tmp(:));
    tmp = r3{w}(tridx{w},2,:);
    ymin(w) = min(tmp(:)); ymax(w) = max(tmp(:));
    tmp = r3{w}(tridx{w},3,:);
    zmin(w) = min(tmp(:)); zmax(w) = max(tmp(:));
    v = zeros(length(tridx{w}),3);
    ww = zeros(length(tridx{w}),3);
    for i = 1:length(tridx{w})
        r2 = projection2(squeeze(r3{w}(tridx{w}(i),:,:)),handles.calib);
        v(i,:) = r2(1,:,2);
        ww(i,:) = r2(2,:,2);
    end
    vmin(w) = min(v(:)); vmax(w) = max(v(:));
    wmin(w) = min(ww(:)); wmax(w) = max(ww(:));
    clear i r2 v ww
end

xmin = min(xmin); xmax = max(xmax);
ymin = min(ymin); ymax = max(ymax);
zmin = min(zmin); zmax = max(zmax);
vmin = min(vmin); vmax = max(vmax);
wmin = min(wmin); wmax = max(wmax);

tridx_all = unique([tridx{:}]);
clear tridx
    
% colours = {'b.-','g.-','r.-','c.-','y.-','m.-'};
% lines = {'b-','g-','r-','c-','y-','m-'};
% points = {'b.','g.','r.','c.','y.','m.'};
t = 0:0.02:1;

for fridx = 1:length(tridx_all)
    fr = tridx_all(fridx);
    clf
    for w = 1:handles.Nwhiskers
        if ~handles.whisker(w).tracked(fr)
            break
        end
        wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
        r3w = squeeze(r3{w}(fr,:,:));
        r2 = projection2(r3w,handles.calib);
        b3 = bezierval(r3w,t);
        b2 = projection2(b3,handles.calib);
        subplot 131,  
        plot3(r3w(1,:),r3w(2,:),r3w(3,:),handles.points{wmod},b3(1,:),b3(2,:),b3(3,:),handles.lines{wmod})
        title(fr), xlabel('x'),ylabel('y'),zlabel('z')
        xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);
        hold on
        set(gca,'ydir','reverse'), 
%         subplot 132, hold on
%         plot(r2(1,:,1),r2(2,:,1),points{wmod},b2(1,:,1),b2(2,:,1),lines{wmod})
%         title('H view'), set(gca,'ydir','reverse'), axis image
%         xlim([xmin xmax]); ylim([ymin ymax]);
%         subplot 133, hold on
%         plot(r2(1,:,2),r2(2,:,2),points{wmod},b2(1,:,2),b2(2,:,2),lines{wmod})
%         title('V view'), set(gca,'ydir','reverse'), axis image
%         xlim([vmin vmax]); ylim([wmin wmax]);
    end
    pause(tpause)
    clear w wmod r3w r2 b3 b2
end
clear fridx fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_number_whiskers_Callback(hObject, eventdata, handles)
% hObject    handle to edit_number_whiskers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_number_whiskers as text
%        str2double(get(hObject,'String')) returns contents of edit_number_whiskers as a double

handles.Nwhiskers = str2double(get(hObject,'String'));
handles.selected_whiskers = ones(1:handles.Nwhiskers);
handles.current_whisker = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_number_whiskers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_number_whiskers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_current_whisker_selected.
function radiobutton_current_whisker_selected_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_current_whisker_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_current_whisker_selected

handles.whisker(handles.current_whisker).selected = get(hObject,'Value');
guidata(hObject, handles);

function edit_current_whisker_label_Callback(hObject, eventdata, handles)
% hObject    handle to edit_current_whisker_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_current_whisker_label as text
%        str2double(get(hObject,'String')) returns contents of edit_current_whisker_label as a double

handles.whisker(handles.current_whisker).label = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_current_whisker_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_current_whisker_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_add_whisker.
function pushbutton_add_whisker_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_whisker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.Nwhiskers = handles.Nwhiskers + 1;
handles.current_whisker = handles.Nwhiskers;
set(handles.text_Nwhiskers,'String',handles.Nwhiskers);
handles.whisker(handles.current_whisker).label = 'xx';
handles.whisker(handles.current_whisker).energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;

% handles.whisker(handles.current_whisker).energy_threshold = 70;
handles.whisker(handles.current_whisker).selected = 1;
update_current_whisker_display(handles.current_whisker, handles);

% initialise
handles.whisker(handles.current_whisker).tracked = zeros(1,handles.video.h.nframes);
handles.whisker(handles.current_whisker).r3all = zeros(handles.video.h.nframes,3,handles.n_control_pts);
handles.whisker(handles.current_whisker).fp3_all = zeros(handles.video.h.nframes,3);
handles.whisker(handles.current_whisker).kappa_all = zeros(handles.video.h.nframes,2);
handles.whisker(handles.current_whisker).theta_all = zeros(handles.video.h.nframes,2);
handles.whisker(handles.current_whisker).Emin = zeros(handles.video.h.nframes,1);
handles.whisker(handles.current_whisker).fpidx = zeros(handles.video.h.nframes,1);
% handles.whisker(handles.current_whisker).folliclemask = cell(handles.video.h.nframes,2);

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_current_whisker_display(w, handles)

set(handles.text_current_whisker,'String',w)
set(handles.edit_current_whisker_label,'String',handles.whisker(w).label)
set(handles.edit_current_whisker_energy_threshold,'String',handles.whisker(w).energy_threshold)
set(handles.radiobutton_current_whisker_selected,'value',handles.whisker(w).selected)

% --- Executes on selection change in popupmenu_set_roi_from_file.
function popupmenu_set_roi_from_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_set_roi_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_set_roi_from_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_set_roi_from_file

val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    trfname = string_list{val};   % avi file or dat file.
else
    error('')
end
clear val string_list

roi_old = handles.roi;
load(trfname,'roi','-mat')
clear trfname

if exist('roi','var')
    save(handles.trfname,'roi','-append')
    handles.roi = roi;
    clear roi
else
    error('File does not contain ''roi'' ')
end

% pmtrs that depend on roi need to be reset:
% handles.video.width.h = handles.roi.h(2)-handles.roi.h(1)+1;
% handles.video.width.v = handles.roi.v(2)-handles.roi.v(1)+1;
% handles.video.height.h = handles.roi.h(4)-handles.roi.h(3)+1;
% handles.video.height.v = handles.roi.v(4)-handles.roi.v(3)+1;

handles.meanframe = init_subtract_image_mean(get(handles.radiobutton_subtract_image_mean,'value'),handles.roi,handles.video);

% load and plot the first frame:
titles.h = ''; titles.v = '';
handles.frame = load_and_plot_frame(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,handles);
clear titles

handles.whisker = translate_contours_to_new_roi(handles.whisker,handles.roi,roi_old);
clear roi_old

% Replot the tracking solution (if any)
if isfield(handles, 'calib')
    plot_contours(handles.currentframe,handles)
end


guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function popupmenu_set_roi_from_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_set_roi_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.tr4');
clear ffavi ffdat

nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Set ROIs from .tr4 file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = find_contour_wm4(r3,im)
% 290116 - this function needs logical rethink.
debug = 0;

D = size(r3,2); % number of control points

% midpoint on P0 - PD-1:
mp = .5*r3(:,1)+.5*r3(:,D);

% normal to P0 - PD-1 at MP:
R = [0 -1;1 0]; % 90deg rotation matrix
nvec = R*(r3(:,D)-rold(:,1));
nvec = nvec/norm(nvec);

% image intensity along nvec
range = 50;
z = -range:range;
v = mp*ones(size(z)) + nvec*z;
i = interp3(im.s,v(1,:),v(2,:),v(3,:));
! this can't work.  this function needs rethinking

% find local minima
di = [0 diff(i)];
minzidx = find((di(1:end-1)<0)&(di(2:end)>0));
if isempty(minzidx)
    error('This should not happen')
end
% find deepish local minima:
% assume that more background than whisker along i, so that median(i) ~
% background
theta = (median(i)+min(i))/2;
minzidx = minzidx(i(minzidx) < theta);
if isempty(minzidx)
    error
end
% find local minimum nearest middle of range:
[~,idx] = min(abs(z(minzidx)));
zstar = z(minzidx(idx));
clear idx

% finally, construct r from zstar
r = zeros(2,3);
r(:,2) = mp + zstar*nvec;
r(:,3) = r(:,2) + .5*(rold(:,3)-rold(:,1));
r(:,1) = mp - .5*(rold(:,3)-rold(:,1));

if debug
    figure
    subplot 211
    imagesc(im.s), colormap gray
    hold on
    plot(rold(1,:),rold(2,:),'r.',mp(1),mp(2),'ro')
    plot(v(1,:),v(2,:),'b')
    plot(r(1,:),r(2,:),'mo')
    subplot 212
    plot(z,i,'k',z(minzidx),i(minzidx),'r.',zstar,50,'r*')
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = z_from_xyvw(xy,vw, M, O)
% Recover 'z' coord from x,y,v,w and calibration
% let P = (x,y,z)'     ie real world coords of a point in 3D
% xy = (x,y)' is projection of P in the horizontal image plane
% vw = (v,w)' is projection of P in the ~vertical image plane
% (NB that xy,vw are projections of the *same* point)
% The calibration equations are:
% v = Mv'*P + 0v
% w = Mw'*P + Ow
% v ~ y, it is w that contains most z information.  So, to get z, invert
% the w equation:

Mw = M(2,:);
Ow = O(2);
w = vw(2);
z = -(Mw(1:2)*xy + Ow - w)/Mw(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_calibrate.
function pushbutton_calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


calib = calibrate_gui_wm4b1(handles.roi, handles.fname, handles.frame, handles.dir);
save(handles.trfname,'calib','-append')
% keyboard

% debugging hack:
% load(handles.trfname,'calib','-mat')

handles.calib = calib;
clear calib

% 060416 is following necessary?
% re-initialise video: 160816 i've commented it out (following lin)
% handles = initialise_new_video(handles.fname,handles);
% handles.video.width.h = handles.roi.h(2)-handles.roi.h(1)+1;
% handles.video.width.v = handles.roi.v(2)-handles.roi.v(1)+1;
% handles.video.height.h = handles.roi.h(4)-handles.roi.h(3)+1;
% handles.video.height.v = handles.roi.v(4)-handles.roi.v(3)+1;

% load and plot the first frame:
titles.h = ''; titles.v = '';
haxes.h = handles.viewh; 
haxes.v = handles.viewv; 
handles.frame = load_and_plot_frames(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,haxes);

% plot contours (using new calib):
plot_contours(handles.currentframe,handles)
clear haxes

set(handles.pushbutton_fit_snakes,'enable','on')

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_calibrate_from_file.
function popupmenu_calibrate_from_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_calibrate_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_calibrate_from_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_calibrate_from_file


val = get(hObject,'Value');
if(val~=1)
    string_list = get(hObject,'String');
    trfname = string_list{val};   % avi file or dat file.
else
    error('')
end
clear val string_list

roi_old = handles.roi;
load(trfname,'calib','-mat')
clear trfname

handles.calib = calib;

% load and plot the first frame:
% handles.frame = load_and_plot_frame(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,handles);
handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);

cla(handles.viewh,'reset')
cla(handles.viewv,'reset')

titles.h = num2str(handles.currentframe.h); titles.v = num2str(handles.currentframe.v);
haxes.h = handles.viewh;
haxes.v = handles.viewv;

plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
clear titles haxes

if isfield(handles,'calib')
%     plot_contours(handles.currentframe.h,handles)
    set(handles.pushbutton_fit_snakes,'enable','on')
end

save(handles.trfname,'-mat','-append','calib')
clear calib

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function popupmenu_calibrate_from_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_calibrate_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.tr4');

nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Set calibration from .tr4 file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pv,pw,z] = hview_point_2_vview_line(x,y,calib,roi)

% The aim here is to take a point (x,y) in H view and calculate the
% corresponding line in V view.  For theory of this, see notebook entry 110216
% The output is polyval coefficients of a line, parameterised as a function
% of 'z'
% ie v(z) = a*z + b, w(z) = c*z + d
% for z in a range defined to match the size of the V view image

% coefficents of linear polynomial in V view:
pv = zeros(1,2);
pv(1) = calib.mv(3);
pv(2) = calib.mv(1)*x + calib.mv(2)*y + calib.ov;
pw = zeros(1,2);
pw(1) = calib.mw(3);
pw(2) = calib.mw(1)*x + calib.mw(2)*y + calib.ow;

zmax = (1 - pw(2))/pw(1);
zmin = (roi.v(4) - pw(2))/pw(1);
z = linspace(zmin,zmax,100);
clear zmin zmax



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma_prior as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma_prior as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma_prior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function edit_energy_threshold_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_energy_threshold (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_energy_threshold as text
% %        str2double(get(hObject,'String')) returns contents of edit_energy_threshold as a double
% 
% handles.energy_threshold = str2double(get(hObject,'String'));
% guidata(hObject, handles);

% % --- Executes during object creation, after setting all properties.
% function edit_energy_threshold_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_energy_threshold (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 


function edit_playback_pause_Callback(hObject, eventdata, handles)
% hObject    handle to edit_playback_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_playback_pause as text
%        str2double(get(hObject,'String')) returns contents of edit_playback_pause as a double


% --- Executes during object creation, after setting all properties.
function edit_playback_pause_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_playback_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0.1')
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_current_whisker_energy_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_current_whisker_energy_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_current_whisker_energy_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_current_whisker_energy_threshold as a double

handles.whisker(handles.current_whisker).energy_threshold = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_current_whisker_energy_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_current_whisker_energy_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_direction_video.
function pushbutton_direction_video_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_direction_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'String')
    case 'forward'
        set(hObject,'String','reverse')
    case 'reverse'
        set(hObject,'String','forward')
    otherwise
        error('Should not happen')
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function pushbutton_direction_video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_direction_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'String','reverse')
guidata(hObject, handles)



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit_playback_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_playback_pause as text
%        str2double(get(hObject,'String')) returns contents of edit_playback_pause as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_playback_pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_start_video.
function radiobutton_start_video_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_start_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_start_video


set(handles.radiobutton_stop_video,'value',false)
play = true;
titles.h = ''; titles.v = '';
dt = 0.05;
t = 0:dt:1;

tmp2.x = []; tmp2.y = []; tmp2.z = [];
for w = 1:handles.Nwhiskers
    tmp.x = handles.whisker(w).r3all(:,1,:);
    tmp.y = handles.whisker(w).r3all(:,2,:);
    tmp.z = handles.whisker(w).r3all(:,3,:);
    tmp2.x = [tmp2.x; tmp.x(:)]; tmp2.y = [tmp2.y; tmp.y(:)]; tmp2.z = [tmp2.z; tmp.z(:)];
end
xmin = min(tmp2.x); xmax = max(tmp2.x);
ymin = min(tmp2.y); ymax = max(tmp2.y);
zmin = min(tmp2.z); zmax = max(tmp2.z);
clear tmp tmp2

while(play)
    tic,
    % Check if user has pressed any relevant buttons:
    if get(handles.radiobutton_stop_video,'value')
        play = false;
    end
    tpause = str2double(get(handles.edit_playback_pause,'String'));
    finterval = str2double(get(handles.edit_videoplayback_frame_interval,'String'));
    switch get(handles.pushbutton_direction_video,'String')
        case 'reverse'
            handles.currentframe.h = modadd(handles.currentframe.h,finterval,handles.video.h.nframes);
            handles.currentframe.v = modadd(handles.currentframe.v,finterval,handles.video.v.nframes);
        case 'forward'
            handles.currentframe.h = modsubtract(handles.currentframe.h,finterval,handles.video.h.nframes);
            handles.currentframe.v = modsubtract(handles.currentframe.v,finterval,handles.video.v.nframes);
        otherwise
            error('Should not happen')
    end
    % Load and plot current video frame:
%     handles.frame = load_frame_2views(handles.video,handles.currentframe,handles.roi);
        handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
        handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);

%     set(handles.currentframe_display,'string',num2str(handles.currentframe.h))
    cla(handles.viewh,'reset')
%     cla(handles.viewv,'reset')
    titles.h = num2str(handles.currentframe.h); titles.v = num2str(handles.currentframe.v);
    haxes.h = handles.viewh;
    haxes.v = handles.viewv;
%     if get(handles.radiobutton_show_3d,'value')
% %         cla(handles.axes_misc_plot,'reset')
%         cla(handles.axes_misc_plot)
%     end
%     plot_frame(handles.frame,handles.video,handles.meanframe,titles,handles)
%             plot_frame(handles.frame,handles.video,handles.meanframe,titles,handles)
    plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
    clear titles haxes

    % Superimpose tracking solutions, if any:
    for w = 1:handles.Nwhiskers
        wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
        if any(handles.whisker(w).tracked)
            r3 = squeeze(handles.whisker(w).r3all(handles.currentframe.h,:,:));
            b3 = bezierval(r3,t);
            r2 = projection2(r3,handles.calib);
            b2 = projection2(b3,handles.calib);
            axes(handles.viewh), hold on
            plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod})
            if ~handles.tracking2D
            axes(handles.viewv), hold on
            plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod})
            end
        end
%         if get(handles.radiobutton_show_3d,'value')
%             axes(handles.axes_misc_plot),
%             plot3(b3(1,:),b3(2,:),b3(3,:),handles.lines{wmod})
%             xlim([xmin xmax])
%             ylim([ymin ymax])
%             zlim([zmin zmax])
%             hold on
%             set(gca,'ydir','reverse')
%         end
    end
    clear w r3 b3 r2 br2 wmod

    pause(tpause)
    
end
guidata(hObject, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in radiobutton_stop_video.
function radiobutton_stop_video_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_stop_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_stop_video

set(handles.radiobutton_start_video,'value',false)
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function radiobutton_start_video_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_start_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in radiobutton_show_3d.
function radiobutton_show_3d_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_show_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_show_3d


% --- Executes on selection change in popupmenu_initialise_from_file.
function popupmenu_initialise_from_file_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_initialise_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_initialise_from_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_initialise_from_file

contents = cellstr(get(hObject,'String'));
trfile = contents{get(hObject,'Value')};
clear contents
load(trfile,'whisker','roi','calib','trpmtrs','-mat')

% Set ROI using the tr file:
if exist('roi','var')
    handles.roi = roi;
    clear roi
    % pmtrs that depend on roi need to be reset:
%     handles.video.width.h = handles.roi.h(2)-handles.roi.h(1)+1;
%     handles.video.width.v = handles.roi.v(2)-handles.roi.v(1)+1;
%     handles.video.height.h = handles.roi.h(4)-handles.roi.h(3)+1;
%     handles.video.height.v = handles.roi.v(4)-handles.roi.v(3)+1;
    handles.meanframe = init_subtract_image_mean(false,handles.roi,handles.video);
%     handles.whisker = translate_contours_to_new_roi(handles.whisker,handles.roi,roi_old);
%     clear roi_old
end

% Set Calib using the tr file:
if exist('calib','var')
    handles.calib = calib;
    clear calib
    set(handles.pushbutton_fit_snakes,'enable','on') 
end

% load and plot the first frame:
% titles.h = ''; titles.v = '';
% handles.frame = load_and_plot_frame(handles.video,handles.currentframe,handles.roi,handles.meanframe,titles,handles);
% clear titles

handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
cla(handles.viewh,'reset')
titles.h = num2str(handles.currentframe.h); 
haxes.h = handles.viewh;

if ~handles.tracking2D
handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);
cla(handles.viewv,'reset')
titles.v = num2str(handles.currentframe.v);
haxes.v = handles.viewv;
end
plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
clear titles haxes


% Initialise whisker parameters from the tr file:
if exist('whisker','var')
    handles.Nwhiskers = length(whisker);
    for w = 1:handles.Nwhiskers
        handles.whisker(w).label = whisker(w).label;
        handles.whisker(w).energy_threshold = whisker(w).energy_threshold;
        handles.whisker(w).selected = whisker(w).selected;
    end
    set(handles.text_Nwhiskers,'String',handles.Nwhiskers)
    handles.current_whisker = 1;
    update_current_whisker_display(handles.current_whisker, handles);
else
    % parameters will be those set up by prior call to initialise_new_track
end
clear trfile whisker
% initialise tracking parameters from the  tr file:
if exist('trpmtrs','var')
    if ~isfield(trpmtrs,'tracking_direction')   % new in v4b2 050916
        % for backwards compatibility with files that might lack this pmtr
        trpmtrs.tracking_direction = handles.trpmtrs.tracking_direction;
    end
    handles.trpmtrs = trpmtrs;
    clear trpmtrs
    set(handles.edit_sigma_prior,'String',handles.trpmtrs.sigma_prior);
    if isfield(handles,'trpmtrs_sigma2_prior')
        set(handles.edit_sigma2_prior,'String',handles.trpmtrs_sigma2_prior);
    end
    if isfield(handles,'trpmtrs_snout_sigma')
        set(handles.edit_snout_sigma,'String',handles.trpmtrs_snout_sigma);
    end
    %set(handles.radiobutton_subtract_image_mean,'value',handles.trpmtrs.subtract_image_mean)
    %set(handles.pushbutton_tracking_direction,'String',handles.trpmtrs_tracking_direction);
end

% initialise tracking variables:
handles.n_control_pts = 3;
for w = 1:handles.Nwhiskers
    handles.whisker(w).tracked = zeros(1,handles.video.h.nframes);
    handles.whisker(w).r3all = zeros(handles.video.h.nframes,3,handles.n_control_pts);
    handles.whisker(w).fp3_all = zeros(handles.video.h.nframes,3);
%     handles.whisker(w).kappa_all = zeros(handles.video.nframes,2);
%     handles.whisker(w).theta_all = zeros(handles.video.nframes,2);
    handles.whisker(w).Emin = zeros(handles.video.h.nframes,1);
    handles.whisker(w).fpidx = zeros(handles.video.h.nframes,1);
    handles.whisker(w).s0 = [];
end
clear w
%set(handles.pushbutton_start_batch,'enable','on')

guidata(hObject, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_tracks_from_file(trfile,handles)
% superimpose solutions from tr file on image plots to guide initialisation:

load(trfile,'whisker','calib','-mat')
if ~exist('whisker','var') | ~exist('calib','var')
    return
end
Nwhiskers = length(whisker);

for w = 1:Nwhiskers
    wmod = rem(w,Nwhiskers) + Nwhiskers*(w==Nwhiskers); % to avoid running out of colours
    trframes = find(whisker(w).tracked);
    if numel(trframes)>20
        trframes = round(linspace(trframes(1),trframes(end),20));
    end
    for fr = trframes
        r3 = squeeze(whisker(w).r3all(fr,:,:));
        b3 = bezierval(r3,0:.02:1);
        b2 = projection2(b3,calib);
        axes(handles.viewh)
        plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod}, 'LineWidth',2)
        if ~handles.tracking2D
        axes(handles.viewv)
        plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod}, 'LineWidth',2),
        end
    end
end
clear w wmod trframes fr r3 b3 b2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function popupmenu_initialise_from_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_initialise_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ff = dir('*.tr4');
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Initialisation file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(hObject,'String',string_list);
clear ff nr_files string_list i


% --- Executes on button press in pushbutton_reset_current_whisker.
function pushbutton_reset_current_whisker_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset_current_whisker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.whisker(handles.current_whisker).tracked = zeros(1,handles.video.h.nframes);
handles.whisker(handles.current_whisker).r3all = zeros(handles.video.h.nframes,3,handles.n_control_pts);
handles.whisker(handles.current_whisker).fp3_all = zeros(handles.video.h.nframes,3);
handles.whisker(handles.current_whisker).s0 = [];
handles.whisker(handles.current_whisker).selected = 1;

whisker = handles.whisker; %#ok<NASGU>
save(handles.trfname,'whisker','-append');

frameidx = handles.currentframe.h;
titles.h = ''; titles.v = '';
% plot_frame(handles.frame,handles.video,handles.meanframe,titles,handles)
haxes.h = handles.viewh;
haxes.v = handles.viewv;
plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
for w = 1:handles.Nwhiskers
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx)
        %         hold on
        %         rall = handles.rall.h;
        dt = 0.05;
        t = 0:dt:1;
        r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
        b3 = bezierval(r3,t);
        r2 = projection2(r3,handles.calib);
        b2 = projection2(b3,handles.calib);
        fp3 = handles.whisker(w).fp3_all(frameidx,:)';
        fp2 = projection2(fp3,handles.calib);
        axes(handles.viewh), hold on
        plot(b2(1,:,1),b2(2,:,1),handles.lines{w},r2(1,:,1),r2(2,:,1),handles.points{w},...
            fp2(1,1),fp2(2,1),'y.')
        if ~handles.tracking2D
            axes(handles.viewv), hold on
            plot(b2(1,:,2),b2(2,:,2),handles.lines{w},r2(1,:,2),r2(2,:,2),handles.points{w},...
            fp2(1,2),fp2(2,2),'y.')
        end
    end
end


guidata(hObject, handles);


% --- Executes on button press in pushbutton_delete_current_whisker.
function pushbutton_delete_current_whisker_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete_current_whisker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Nwhiskers>1
    whiskers_to_keep = setdiff(1:handles.Nwhiskers,handles.current_whisker);
    handles.Nwhiskers = numel(whiskers_to_keep);
    handles.whisker = handles.whisker(1,whiskers_to_keep);
    handles.current_whisker = 1;
else
    % handles.Nwhiskers==1, current_whisker==1
    % reset
    handles.whisker(1).tracked = zeros(1,handles.video.h.nframes);
    handles.whisker(1).r3all = zeros(handles.video.h.nframes,3,handles.n_control_pts);
    handles.whisker(1).fp3_all = zeros(handles.video.h.nframes,3);
    handles.whisker(1).s0 = [];
    handles.whisker(1).selected = 1;
    handles.whisker(1).label = 'xx';
    if get(handles.radiobutton_subtract_image_mean,'Value')
        handles.whisker(1).energy_threshold = handles.default_energy_threshold_subtract_image_mean;
    else
        handles.whisker(1).energy_threshold = handles.default_energy_threshold_do_not_subtract_image_mean;
    end
end
% update gui display:
set(handles.text_Nwhiskers,'String',handles.Nwhiskers)
set(handles.text_current_whisker,'String',handles.current_whisker)
set(handles.edit_current_whisker_label,'String',handles.whisker(handles.current_whisker).label)
set(handles.radiobutton_current_whisker_selected,'value',handles.whisker(handles.current_whisker).selected)
set(handles.edit_current_whisker_energy_threshold,'String', handles.whisker(handles.current_whisker).energy_threshold);

whisker = handles.whisker;
save(handles.trfname,'whisker','-append');
clear whisker

frameidx = handles.currentframe;
titles.h = ''; titles.v = '';
haxes.h = handles.viewh;
haxes.v = handles.viewv;
plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)

for w = 1:handles.Nwhiskers
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx)
        %         hold on
        %         rall = handles.rall.h;
        dt = 0.05;
        t = 0:dt:1;
        r3 = squeeze(handles.whisker(w).r3all(frameidx,:,:));
        b3 = bezierval(r3,t);
        r2 = projection2(r3,handles.calib);
        b2 = projection2(b3,handles.calib);
        fp3 = handles.whisker(w).fp3_all(frameidx,:)';
        fp2 = projection2(fp3,handles.calib);
        axes(handles.viewh), hold on
        plot(b2(1,:,1),b2(2,:,1),handles.lines{w}, 'LineWidth',2,r2(1,:,1),r2(2,:,1),handles.points{w},...
            fp2(1,1),fp2(2,1),'y.')
        if ~handles.tracking2D
            axes(handles.viewv), hold on
            plot(b2(1,:,2),b2(2,:,2),handles.lines{w}, 'LineWidth',2,r2(1,:,2),r2(2,:,2),handles.points{w},...
            fp2(1,2),fp2(2,2),'y.')
        end
    end
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton_select_all_whiskers.
function pushbutton_select_all_whiskers_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select_all_whiskers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for w = 1:handles.Nwhiskers
    handles.whisker(w).selected = 1;
end
clear w
set(handles.radiobutton_current_whisker_selected,'value',handles.whisker(handles.current_whisker).selected)

frameidx = handles.currentframe;
titles.h = ''; titles.v = '';
haxes.h = handles.viewh;
haxes.v = handles.viewv;
plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
for w = 1:handles.Nwhiskers
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours
    if handles.whisker(w).selected && handles.whisker(w).tracked(frameidx.h)
        dt = 0.05;
        t = 0:dt:1;
        r3 = squeeze(handles.whisker(w).r3all(frameidx.h,:,:));
        b3 = bezierval(r3,t);
        r2 = projection2(r3,handles.calib);
        b2 = projection2(b3,handles.calib);
        fp3 = handles.whisker(w).fp3_all(frameidx.h,:)';
        fp2 = projection2(fp3,handles.calib);
        axes(handles.viewh), hold on
        %plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod},)
        plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod},r2(1,:,1),r2(2,:,1),handles.points{wmod},...
            fp2(1,1),fp2(2,1),'y.')
        if ~handles.tracking2D
            axes(handles.viewv), hold on
            plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod},r2(1,:,2),r2(2,:,2),handles.points{wmod},...
            fp2(1,2),fp2(2,2),'y.')
        end
    end
end

guidata(hObject, handles);


% --- Executes on button press in radiobutton_plot_snout_contour.
function radiobutton_plot_snout_contour_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_plot_snout_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_plot_snout_contour



function edit_videoplayback_frame_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_videoplayback_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_videoplayback_frame_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_videoplayback_frame_interval as a double


% --- Executes during object creation, after setting all properties.
function edit_videoplayback_frame_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_videoplayback_frame_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_goto_start.
function pushbutton_goto_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_goto_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.radiobutton_stop_video,'value','true')

handles.currentframe.h = handles.video.h.startframe;
handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
titles.h = num2str(handles.currentframe.h);
haxes.h = handles.viewh;

if ~handles.tracking2D
    handles.currentframe.v = handles.video.v.startframe;
    handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);
    titles.v = num2str(handles.currentframe.v);
    haxes.v = handles.viewv;

end

plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
clear titles haxes

guidata(hObject, handles)


% --- Executes on button press in pushbutton_tracking_direction.
function pushbutton_tracking_direction_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tracking_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


switch handles.trpmtrs.tracking_direction
    case 'fwd'
        handles.trpmtrs.tracking_direction = 'bkwd';
    case 'bkwd'
        handles.trpmtrs.tracking_direction = 'fwd';
    otherwise
        error(['Unrecognised tracking direction ' handles.trpmtrs.tracking_direction])
end

set(hObject,'String',handles.trpmtrs.tracking_direction)
guidata(hObject, handles)



function edit_snout_outliers_Callback(hObject, eventdata, handles)
% hObject    handle to edit_snout_outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_snout_outliers as text
%        str2double(get(hObject,'String')) returns contents of edit_snout_outliers as a double


% --- Executes during object creation, after setting all properties.
function edit_snout_outliers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_snout_outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



    function plot_kinematics_test(handles,w,wmod,ttest)
        % test curvature
        % plot circle with given curvature
        r3 = squeeze(handles.whisker(w).r3all(handles.currentframe.h,:,:));
        r2 = projection2(r3,handles.calib);
        % compute curvature projected onto image planes (what we need here)
        kappa.h = curvature(squeeze(r2(:,:,1)),ttest);
        kappa.v = curvature(squeeze(r2(:,:,2)),ttest);
        %             clear r
        rad.h = abs(1/kappa.h);
        rad.v = abs(1/kappa.v);
        %             r = squeeze(handles.whisker(w).r3all(handles.currentframe.h,:,:));
        cc3 = bezierval(r3,ttest);
        cc2 = projection2(cc3,handles.calib);
        cc.h = cc2(:,1);
        cc.v = cc2(:,2);
        clear cc3
        tmp3 = bezierval(r3,ttest+[-.01 .01]);
        tmp2 = projection2(tmp3,handles.calib);
        clear r ttest
        %             tv.h = diff(tmp([1 2],:),1,2); tv.h = tv.h/norm(tv.h);
        %             tv.v = diff(tmp([2 3],:),1,2); tv.v = tv.v/norm(tv.v);
        tv.h = diff(squeeze(tmp2(:,:,1)),1,2); tv.h = tv.h/norm(tv.h);
        tv.v = diff(squeeze(tmp2(:,:,2)),1,2); tv.v = tv.v/norm(tv.v);
        clear tmp
        %             keyboard
        nv.h = [0 1;-1 0] * tv.h;
        nv.v = [0 1;-1 0] * tv.v;
        clear tv
        cc.h = cc.h + sign(-kappa.h) * rad.h * nv.h;
        cc.v = cc.v + sign(-kappa.v) * rad.v * nv.v;
        clear kappa nv
        tall = 0:.01:2*pi;
        c.h = zeros(2,length(tall));
        c.v = zeros(2,length(tall));
        for tidx = 1:length(tall)
            t = tall(tidx);
            %                 rad*[cos(t);sin(t)],cc
            c.h(:,tidx) = cc.h + rad.h*[cos(t);sin(t)];
            c.v(:,tidx) = cc.v + rad.v*[cos(t);sin(t)];
        end
        clear rad cc tall tidx t
        axes(handles.viewh)
        plot(c.h(1,:),c.h(2,:),handles.lines_dotted{wmod})
        if ~handles.tracking2D
        axes(handles.viewv)
        plot(c.v(1,:),c.v(2,:),handles.lines_dotted{wmod})
        end
        clear c



function edit_kinematics_test_position_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kinematics_test_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kinematics_test_position as text
%        str2double(get(hObject,'String')) returns contents of edit_kinematics_test_position as a double


% --- Executes during object creation, after setting all properties.
function edit_kinematics_test_position_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kinematics_test_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_tracking_2D.
function pushbutton_tracking_2D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tracking_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disable vertical view and related buttons
set(handles.pushbutton_calibrate,'enable','off')
set(handles.popupmenu_choose_vvideo,'enable','off')
set(handles.popupmenu_calibrate_from_file,'enable','off')
set(handles.viewv,'Visible','off')
handles.tracking2D=true;

d = handles.dir.h;
ff = dir([d '*.tr4_2D']);
clear d
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Initialisation file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_initialise_from_file,'String',string_list);
set(handles.popupmenu_choose_auto_initialise_file,'String',string_list);
clear nr_files string_list i

guidata(hObject, handles);


% --- Executes on button press in pushbutton_tracking_3D.
function pushbutton_tracking_3D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tracking_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_calibrate,'enable','on')
set(handles.popupmenu_choose_vvideo,'enable','on')
set(handles.popupmenu_calibrate_from_file,'enable','on')
set(handles.viewv,'Visible','on')
handles.tracking2D=false;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_replicate_last_f.
function pushbutton_replicate_last_f_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_replicate_last_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%copy all  the selected whiskers
switch handles.trpmtrs.tracking_direction   % v42b 050916
    case 'fwd'
        track_fwd = true;
    case 'bkwd'
        track_fwd = false;
    otherwise
        error(['Unrecognised handles.trpmtrs.tracking_direction ' handles.trpmtrs.tracking_direction])
end

D = handles.n_control_pts;
% h video is master, v video slave.  So, eg, 'lastframe' is lastframe of h video.
if track_fwd   % v42b 050916
    lastframe = modsubtract(handles.currentframe.h,1,handles.video.h.nframes);
    lastbutoneframe = modsubtract(handles.currentframe.h,2,handles.video.h.nframes);
else
    lastframe = modadd(handles.currentframe.h,1,handles.video.h.nframes);
    lastbutoneframe = modadd(handles.currentframe.h,2,handles.video.h.nframes);
end


%%check that #is not the first or last frame (no initialization)
for w = 1:handles.Nwhiskers

    if handles.whisker(w).selected
        if (handles.currentframe.h==handles.video.h.startframe) || ~handles.whisker(w).tracked(lastframe)
            disp('First frame (no initialization)')
            return;
        end
    end
end

doplot=1;
doplot_full=1;
doplot_light=0;
for w = 1:handles.Nwhiskers
    if handles.whisker(w).selected==false
        continue
    end
    
    wmod = rem(w,handles.Nwhiskers) + handles.Nwhiskers*(w==handles.Nwhiskers); % to avoid running out of colours

    %%copy snout_contour

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%get follicle mask every 2 frame
    if ~isfield(handles,'folliclemask')
        roi.h.x = 50:size(handles.frame.h,2);%50
        roi.h.y = 50:size(handles.frame.h,1);%50
        [folliclemask{1}] = snout_segment(double(handles.frame.h(:,:,1)),roi.h,[],[],[],'horizontal',handles);
        handles.folliclemask{handles.currentframe.h,1} = folliclemask{1};
        if ~handles.tracking2D
            roi.v.x = 50:size(handles.frame.v,2);%50
            roi.v.y = 50:size(handles.frame.v,1);%50
            [folliclemask{2}] = snout_segment(double(handles.frame.v(:,:,1)),roi.v,[],[],[],'vertical',handles);
            handles.folliclemask{handles.currentframe.h,2} = folliclemask{2};
        end
        clear folliclemask_old fpidx_old roi %theta
    else
        
        [folliclemask{1}]=handles.folliclemask{lastframe,1};
        handles.folliclemask{handles.currentframe.h,1} = folliclemask{1};
        if ~handles.tracking2D
            [folliclemask{2}]=handles.folliclemask{lastframe,2};
            handles.folliclemask{handles.currentframe.h,2} = folliclemask{2};
        end
        
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(handles.radiobutton_plot_snout_contour,'value')
        axes(handles.viewh)
        plot(folliclemask{1}(1,:),folliclemask{1}(2,:),'y-')
        if ~handles.tracking2D
            axes(handles.viewv)
            plot(folliclemask{2}(1,:),folliclemask{2}(2,:),'y-')
        end
    end
    update_current_whisker_display(w, handles)
    r3all = handles.whisker(w).r3all;
    r3new = reshape(r3all(lastframe,:,:),3,D);
    r3=r3new;
    
  %%update the structures and frames
    handles.whisker(w).r3all(handles.currentframe.h,:,:) =  handles.whisker(w).r3all(lastframe,:,:);
    handles.whisker(w).fp3_all(handles.currentframe.h,:,:) = handles.whisker(w).fp3_all(lastframe,:,:);
    handles.whisker(w).tracked(handles.currentframe.h) = true;
    %this is not true, compute the real energy?
    handles.whisker(w).Emin(handles.currentframe.h) = handles.whisker(w).Emin(lastframe);
    Emin=handles.whisker(w).Emin(handles.currentframe.h);
    b3 = bezierval(r3new,0:.01:1);
   
    
    if doplot==doplot_full
        fp3=handles.whisker(w).fp3_all(handles.currentframe.h,:)';
        t = [0:handles.dt:1];
        b3 = bezierval(r3new,t);
        r2new = projection2(r3new,handles.calib);
        fp2 = projection2(fp3,handles.calib);
        b2 = projection2(b3,handles.calib);
        titlestring = sprintf('h-frame %d, whisker %d (%.1f)', handles.currentframe.h, w, Emin);
        axes(handles.viewh)
        plot(b2(1,:,1),b2(2,:,1),handles.lines{wmod},fp2(1,1),fp2(2,1),'y.')
        plot(r2new(1,:,1),r2new(2,:,1),handles.points{wmod}),
        if ~handles.tracking2D
            axes(handles.viewv)
            plot(b2(1,:,2),b2(2,:,2),handles.lines{wmod},fp2(1,2),fp2(2,2),'y.')
            plot(r2new(1,:,2),r2new(2,:,2),handles.points{wmod}),
            set(handles.currentframe_display,'String',titlestring);
            
        end
        clear fp3 t 
    end
    if doplot>=doplot_light
        titlestring = sprintf('h-frame %d, whisker %d (%.1f)', handles.currentframe.h, w, Emin);
        set(handles.currentframe_display,'String',titlestring);
        clear titlestring
    end

end
% Proceed to track next frame or stop?
    % v4b2 050916: important change to logic in following lines
%     if handles.currentframe.h ~= handles.video.h.stopframe;   %v4b1
    onwards = false;
    if track_fwd && (handles.currentframe.h~=handles.video.h.stopframe)
        % we're tracking forwards and it's not the final frame, so carry on
        onwards = true;
        handles.currentframe.h = modadd(handles.currentframe.h,1,handles.video.h.nframes);
        if ~handles.tracking2D
        handles.currentframe.v = modadd(handles.currentframe.v,1,handles.video.v.nframes);
        end
    elseif ~track_fwd && (handles.currentframe.h~=handles.video.h.startframe)
        % we're tracking backwards and it's not the first frame, so carry on
        onwards = true;
        handles.currentframe.h = modsubtract(handles.currentframe.h,1,handles.video.h.nframes);
        if ~handles.tracking2D
        handles.currentframe.v = modsubtract(handles.currentframe.v,1,handles.video.v.nframes);
        end
    end

    if onwards
    %         handles.frame = load_frame_2views(handles.video,handles.currentframe,handles.roi);%,...
        handles.frame.h = load_frame(handles.video.h,handles.currentframe.h,handles.roi.h);%,get(handles.checkbox_rotatevview,'value'));
        if ~handles.tracking2D
        handles.frame.v = load_frame(handles.video.v,handles.currentframe.v,handles.roi.v);
        end
%         set(handles.currentframe_display,'string',num2str(handles.currentframe.h))
        if doplot
            pause(.01)
            cla(handles.viewh,'reset')
            titles.h = num2str(handles.currentframe.h); 
            haxes.h = handles.viewh; 
            
            if ~handles.tracking2D
            cla(handles.viewv,'reset')
            titles.v = num2str(handles.currentframe.v);
            haxes.v = handles.viewv; 
            end
            plot_frames(handles.frame,handles.video,handles.meanframe,titles,haxes)
            clear haxes titles
            
        end
    else
        % end of video, stop
        track = 0;
        beep,pause(.5),beep,pause(.5),beep,pause(.5),beep,pause(.5),beep
    end
    clear onwards




%save in the tr4 file
    whisker = handles.whisker;
    trpmtrs = handles.trpmtrs;
    roi = handles.roi;
    calib = handles.calib;
    save(handles.trfname,'whisker','trpmtrs','roi','calib')
    clear whisker trpmtrs roi calib 
%update handles after the for
guidata(hObject, handles);


% --- Executes on button press in pushbutton_tracking_3D.
function tracking_3D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tracking_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_calibrate,'enable','on')
set(handles.popupmenu_choose_vvideo,'enable','on')
set(handles.popupmenu_calibrate_from_file,'enable','on')
set(handles.viewv,'Visible','on')
handles.tracking2D=false;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_tracking_2D.
function tracking_2D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tracking_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_calibrate,'enable','off')
set(handles.popupmenu_choose_vvideo,'enable','off')
set(handles.popupmenu_calibrate_from_file,'enable','off')
set(handles.viewv,'Visible','off')
handles.tracking2D=true;

d = handles.dir.h;
ff = dir([d '*.tr4_2D']);
clear d
nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Initialisation file';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end
set(handles.popupmenu_initialise_from_file,'String',string_list);
set(handles.popupmenu_choose_auto_initialise_file,'String',string_list);
clear nr_files string_list i

guidata(hObject, handles);
