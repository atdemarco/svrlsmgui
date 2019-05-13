function varargout = svrinteract(varargin)
% SVRINTERACT MATLAB code for svrinteract.fig
%      SVRINTERACT, by itself, creates a new SVRINTERACT or raises the existing
%      singleton*.
%
%      H = SVRINTERACT returns the handle to a new SVRINTERACT or the handle to
%      the existing singleton*.
%
%      SVRINTERACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVRINTERACT.M with the given input arguments.
%
%      SVRINTERACT('Property','Value',...) creates a new SVRINTERACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before svrinteract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to svrinteract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help svrinteract

% Last Modified by GUIDE v2.5 19-Mar-2018 12:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @svrinteract_OpeningFcn, ...
                   'gui_OutputFcn',  @svrinteract_OutputFcn, ...
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


% --- Executes just before svrinteract is made visible.
function svrinteract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to svrinteract (see VARARGIN)

% Choose default command line output for svrinteract
handles.output = hObject;
handles.variables = varargin{1};
handles.opts.stand = false;
handles.opts.ep = .1;
handles.opts.box = 30;
handles.opts.ks = 2;
handles.opts.nfolds = 1; % nfolds = 1 is no crossval
handles.opts.zslice = 25;

tmp = 10 .* randn(handles.variables.vo.dim(1:3));

if (numel(handles.variables.vo.dim) > 3) && (handles.variables.vo.dim(4) > 1) % then it's 4D % handles.variables.vo.dim(4) > 1 % then it's 4d...
    handles.ndims = 4;
else
    handles.ndims = 3;
end

if handles.ndims == 4
    handles.I = imagesc(tmp(:,:,20),'parent',handles.axes1,[-50 50]);
else % it's 3d...
    handles.I = imagesc(tmp(:,:,20),'parent',handles.axes1,[-10 10]);
end

handles.montage_zslices = 5:10:(size(tmp,3)-1);

if handles.ndims == 4
    montage(tmp, 'Indices',handles.montage_zslices,'DisplayRange', [-50 50],'parent',handles.axes1);
else
    montage(tmp, 'Indices',handles.montage_zslices,'DisplayRange', [-10 10],'parent',handles.axes1);
end
 colormap jet;
 colorbar(handles.axes1)

realdata=handles.variables.one_score;
handles.real = plot(1:numel(handles.variables.one_score),realdata,'sg-','parent',handles.axes2);
hold on;
set(gca,'ylim',[min(realdata)-.1*min(realdata) max(realdata)+.1*max(realdata)])
hold on;
handles.pred = plot(1:numel(handles.variables.one_score),zeros(1,numel(handles.variables.one_score)),'rx-','parent',handles.axes2);

% Update handles structure
guidata(hObject, handles);

paint_current(hObject,handles)

% UIWAIT makes svrinteract wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function paint_current(hObject,handles)
    set(handles.standardize_toggle,'value',handles.opts.stand)
    set(handles.epsilon_editbox,'string',num2str(handles.opts.ep)); 
    set(handles.cost_editbox,'string',num2str(handles.opts.box)); 
    set(handles.kenelscale_editbox,'string',num2str(handles.opts.ks)); 
    set(handles.zslice,'string',num2str(handles.opts.zslice)); 
    set(handles.crossval_nfolds_editbox,'string',num2str(handles.opts.nfolds)); 
    
    if handles.ndims == 4 % then use linear kernel
    Mdl = fitrsvm(handles.variables.lesion_dat,handles.variables.one_score,'KernelFunction','linear',... % 'rbf', ...
        'KernelScale',handles.opts.ks,'BoxConstraint',handles.opts.box,'Standardize', ...
        handles.opts.stand,'Epsilon',handles.opts.ep);
    else % rbf
    Mdl = fitrsvm(handles.variables.lesion_dat,handles.variables.one_score,'KernelFunction','rbf', ...
        'KernelScale',handles.opts.ks,'BoxConstraint',handles.opts.box,'Standardize', ...
        handles.opts.stand,'Epsilon',handles.opts.ep);
    end
    w = Mdl.Alpha.'*Mdl.SupportVectors;
    handles.variables.beta_scale = 10/max(abs(w));
    
    if handles.ndims == 4
        disp('4D data input...')
        tmp = zeros(handles.variables.vo.dim(1:4)); % THIS IS DIFFERENT FOR 4D DATA
        beta_map = tmp;
        tmp(handles.variables.l_idx) = w'*handles.variables.beta_scale; % return all lesion data to its l_idx indices.
        beta_map(handles.variables.m_idx) = tmp(handles.variables.m_idx); % m_idx -> m_idx
        
        beta_map = sum(beta_map,4); % THIS IS DIFFERENT FOR 4D DATA
        beta_map_bin = beta_map~=0;

        relslices = find(squeeze(sum(squeeze(sum(beta_map_bin,1)),1)));
        handles.montage_zslices = min(relslices):2:max(relslices);
        montage(beta_map,'indices',handles.montage_zslices, 'DisplayRange', [-20 20],'parent',handles.axes1);
    else
        disp('3D data input...')
        tmp = zeros(handles.variables.vo.dim(1:3));
        beta_map = tmp;

        tmp(handles.variables.l_idx) = w'*handles.variables.beta_scale; % return all lesion data to its l_idx indices.
        beta_map(handles.variables.m_idx) = tmp(handles.variables.m_idx); % m_idx -> m_idx

        %set(handles.I,'cdata',beta_map(:,:,handles.opts.zslice));
        %montage(beta_map, 'Indices',handles.montage_zslices,'DisplayRange', [-10 10],'parent',handles.axes1);

        %handles.montage_zslices = 5:10:(size(tmp,3)-1);
        %handles.variables.lesion_dat

        beta_map_bin = beta_map~=0;

        relslices = find(squeeze(sum(squeeze(sum(beta_map_bin,1)),1)));
        handles.montage_zslices = min(relslices):2:max(relslices);
        montage(beta_map,'indices',handles.montage_zslices, 'DisplayRange', [-10 10],'parent',handles.axes1);

    end


    colormap jet;
    colorbar(handles.axes1)

    disp(['Proportion is support vectors = ' num2str(sum(Mdl.IsSupportVector)/numel(Mdl.IsSupportVector))])

    if handles.opts.nfolds == 1 % no crossval
        predicted = Mdl.predict(handles.variables.lesion_dat);
    else
        XVMdl = crossval(Mdl,'KFold',handles.opts.nfolds);
        predicted = XVMdl.kfoldPredict;        
    end


    set(handles.pred,'ydata',predicted)
    hold on;
    set(handles.real,'ydata',handles.variables.one_score)
%     
%     realdata=handles.variables.one_score;
%     handles.real = plot(1:numel(handles.variables.one_score),realdata,'sg-','parent',handles.axes2);
%     hold on;
%     set(handles.axes2,'ylim',[min(realdata)-.1*min(realdata) max(realdata)+.1*max(realdata)])
%    preddata = Mdl.predict(handles.variables.lesion_dat);
%     handles.pred = plot(1:numel(handles.variables.one_score),preddata,'rx-','parent',handles.axes2);
    guidata(hObject, handles);


function varargout = svrinteract_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function epsilon_editbox_Callback(hObject, eventdata, handles)
    handles.opts.ep = str2double(get(hObject,'String'));
    paint_current(hObject,handles)

% --- Executes during object creation, after setting all properties.
function epsilon_editbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cost_editbox_Callback(hObject, eventdata, handles)
    handles.opts.box = str2double(get(hObject,'String'));
    paint_current(hObject,handles)


% --- Executes during object creation, after setting all properties.
function cost_editbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kenelscale_editbox_Callback(hObject, eventdata, handles)
    handles.opts.ks = str2double(get(hObject,'String'));
    paint_current(hObject,handles)


% --- Executes during object creation, after setting all properties.
function kenelscale_editbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in standardize_toggle.
function standardize_toggle_Callback(hObject, eventdata, handles)
    handles.opts.stand = ~handles.opts.stand;
    paint_current(hObject,handles)


% --- Executes on button press in pushbutton2refre.
function pushbutton2refre_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2refre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function zslice_Callback(hObject, eventdata, handles)
    handles.opts.zslice = str2double(get(hObject,'String'));
    paint_current(hObject,handles)



% --- Executes during object creation, after setting all properties.
function zslice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function crossval_nfolds_editbox_Callback(hObject, eventdata, handles)
    handles.opts.nfolds = str2double(get(hObject,'String'));
    paint_current(hObject,handles)

% --- Executes during object creation, after setting all properties.
function crossval_nfolds_editbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
