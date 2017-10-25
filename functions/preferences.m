function varargout = preferences(varargin)
% PREFERENCES MATLAB code for preferences.fig
%      PREFERENCES, by itself, creates a new PREFERENCES or raises the existing
%      singleton*.
%
%      H = PREFERENCES returns the handle to a new PREFERENCES or the handle to
%      the existing singleton*.
%
%      PREFERENCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREFERENCES.M with the given input arguments.
%
%      PREFERENCES('Property','Value',...) creates a new PREFERENCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preferences_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preferences_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preferences

% Last Modified by GUIDE v2.5 17-Apr-2017 14:31:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preferences_OpeningFcn, ...
                   'gui_OutputFcn',  @preferences_OutputFcn, ...
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


% --- Executes just before preferences is made visible.
function preferences_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preferences (see VARARGIN)

% set the editbox values from input fctn
handles.parameters = varargin{1};
handles.details = varargin{2};
set(handles.costeditbox,'string',num2str(handles.parameters.cost))
set(handles.gammaeditbox,'string',num2str(handles.parameters.gamma))

% Set SVR algorithm radio button values.
set(handles.svralgorithm_matlab_radiobutton,'value',~handles.parameters.useLibSVM)
set(handles.svralgorithm_libsvm_radiobutton,'value',handles.parameters.useLibSVM)
if handles.parameters.useLibSVM
    set(handles.gammaeditbox,'enable','on')
    set(handles.costeditbox,'enable','on')
else
    set(handles.gammaeditbox,'enable','off')
    set(handles.costeditbox,'enable','off')
end

if handles.details.stats_toolbox
    set(handles.svralgorithm_matlab_radiobutton,'enable','on')
else
    set(handles.svralgorithm_matlab_radiobutton,'enable','off')
end

if handles.details.libsvm
    set(handles.svralgorithm_libsvm_radiobutton,'enable','on')
else
    set(handles.svralgorithm_libsvm_radiobutton,'enable','off')
end

% Choose default command line output for preferences
handles.output = ''; % default to blank. (orig is hObject)

% center the fig
figpos = get(0,'defaultfigureposition');
oldunits = get(hObject,'units');
set(hObject,'units','pixels');
oldpos = get(hObject,'position');
fig_width = oldpos(3);
fig_height = oldpos(4);

old_gcbf_units = get(gcbf,'units');
set(gcbf,'units','pixels')
gcbfpos = get(gcbf,'position');
set(gcbf,'units',old_gcbf_units)

left_and_top_offset = [ (sum(gcbfpos([1 3]))/2)  - (fig_width/2) ...
                        (sum(gcbfpos([2 4]))/2)  - (fig_height/2) ];
newfigpos = [left_and_top_offset fig_width fig_height];
set(hObject,'position',newfigpos,'units',oldunits)

set(handles.figure1,'WindowStyle','modal'); % make the gui modal
guidata(hObject, handles); % Update handles structure
uiwait(handles.figure1); % UIWAIT makes preferences wait for user response (see UIRESUME)

% --- Outputs from this function are returned to the command line.
function varargout = preferences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);

function gammaeditbox_Callback(hObject, eventdata, handles)
    str = get(gcbo,'string');
    if isempty(str2num(str)) % then it's not a valid number...
        set(gcbo,'string','5');
        warndlg('Input must be numerical [default gamma = 5]');
    else % update the parameter value.
        handles.parameters.gamma = str2num(get(gcbo,'string'));
    end
    guidata(hObject,handles)
    
function gammaeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function costeditbox_Callback(hObject, eventdata, handles)
    str = get(gcbo,'string');
    if isempty(str2num(str)) % then it's not a valid number...
        set(gcbo,'string','30');
        warndlg('Input must be numerical [default cost = 30]');
    else % update the parameter value.
        handles.parameters.cost = str2num(get(gcbo,'string'));
    end
    guidata(hObject,handles)

function costeditbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
    finish_UI(hObject,eventdata,handles)

% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
    finish_UI(hObject,eventdata,handles)
    
function finish_UI(hObject,eventdata,handles)
switch get(gcbo,'tag')
    case 'cancelbutton', handles.output = []; % cues to parent gui that the dialog was canceled.
    case 'okbutton', handles.output = handles.parameters;
end
guidata(hObject,handles)
uiresume(handles.figure1)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting'), uiresume(hObject)
else, delete(hObject)
end

function svralgorithm_matlab_radiobutton_Callback(hObject, eventdata, handles)
    handles.parameters.useLibSVM = 0;
    % Set SVR algorithm radio button values.
    set(handles.svralgorithm_matlab_radiobutton,'value',~handles.parameters.useLibSVM)
    set(handles.svralgorithm_libsvm_radiobutton,'value',handles.parameters.useLibSVM)
    set(handles.gammaeditbox,'enable','off')
    set(handles.costeditbox,'enable','off')
    guidata(hObject,handles)
    

function svralgorithm_libsvm_radiobutton_Callback(hObject, eventdata, handles)
    handles.parameters.useLibSVM = 1;
    % Set SVR algorithm radio button values.
    set(handles.svralgorithm_matlab_radiobutton,'value',~handles.parameters.useLibSVM)
    set(handles.svralgorithm_libsvm_radiobutton,'value',handles.parameters.useLibSVM)
    set(handles.gammaeditbox,'enable','on')
    set(handles.costeditbox,'enable','on')
    guidata(hObject,handles)
