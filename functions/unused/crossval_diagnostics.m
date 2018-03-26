function varargout = crossval_diagnostics(varargin)
% CROSSVAL_DIAGNOSTICS MATLAB code for crossval_diagnostics.fig
%      CROSSVAL_DIAGNOSTICS, by itself, creates a new CROSSVAL_DIAGNOSTICS or raises the existing
%      singleton*.
%
%      H = CROSSVAL_DIAGNOSTICS returns the handle to a new CROSSVAL_DIAGNOSTICS or the handle to
%      the existing singleton*.
%
%      CROSSVAL_DIAGNOSTICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROSSVAL_DIAGNOSTICS.M with the given input arguments.
%
%      CROSSVAL_DIAGNOSTICS('Property','Value',...) creates a new CROSSVAL_DIAGNOSTICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crossval_diagnostics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crossval_diagnostics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crossval_diagnostics

% Last Modified by GUIDE v2.5 21-Feb-2018 11:09:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crossval_diagnostics_OpeningFcn, ...
                   'gui_OutputFcn',  @crossval_diagnostics_OutputFcn, ...
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


% --- Executes just before crossval_diagnostics is made visible.
function crossval_diagnostics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crossval_diagnostics (see VARARGIN)

% Choose default command line output for crossval_diagnostics
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crossval_diagnostics wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crossval_diagnostics_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
