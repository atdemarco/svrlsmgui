function handles = ConfigureSVRLSMGUIOptions(handles)
%     handles.options.lesionvolumecorrection = {'Regress on Behavior','Regress on Lesion','Regress on Both','DTLVC','None'};
%     handles.options.old_hypodirection = {'One-tailed (positive)','One-tailed (negative)','Two-tailed'}; % before 2/7/18
%     handles.options.hypodirection = {'High scores are bad','High scores are good','Two-tailed'};
    handles.options = lsmtb_options;% this should be made obsolete...
    set(handles.lesionvolumecorrectiondropdown,'String',handles.options.lesionvolumecorrection)
    set(handles.hypodirectiondropdown,'String',handles.options.hypodirection)

    % We disable this because it doesn't produce reliable results at the moment...
    set(handles.optimize_menu,'enable','off')
    set(handles.summary_paramoptimization,'checked','off','enable','off')
    
    %% Add generic callbacks
    add_generic_callback_objs = [handles.lesionthresholdeditbox handles.addcovariate handles.removecovariate  ...
        handles.applycovariatestobehaviorcheckbox handles.applycovariatestolesioncheckbox ...
        handles.chooseoutputfolderbutton handles.chooselesionfolderbutton handles.choosescorefilebutton ...
        handles.clusterwisepeditbox handles.cfwer_p_value_editbox handles.cluster_voxelwise_p_editbox handles.npermutationseditbox ...
        handles.cfwer_v_value_editbox handles.permutationtestingcheckbox handles.do_cfwer_checkbox handles.lesionvolumecorrectiondropdown ...
        handles.hypodirectiondropdown handles.scorenamepopupmenu ...
        handles.open_output_folder_button handles.open_lesion_folder_button handles.open_score_file_button];
    
    set(add_generic_callback_objs,'Callback',@(hObject,eventdata)svrlsmgui('UpdateCurrentAnalysis',guidata(hObject),hObject));
    
    generic_menus = [get(handles.search_strategy_options,'children') ; 
                     get(handles.output_summary_menu,'children') ; 
                     get(handles.save_perm_data,'children') ; 
                     handles.lsm_method_parent_menu ; 
                     get(handles.lsm_method_parent_menu,'children') ; 
                     get(handles.search_strategy_menu_option,'children') ;
                     get(handles.parameters_to_optimize_menu,'children') ;
                     get(handles.crossvalidation_parent_menu,'children') ;
                     get(handles.parent_cache_menu,'children') ;
                     handles.optimization_is_verbose_menu ;
                     get(handles.beta_options_menu,'children');
                     handles.image_data_options_parent_menu;
                     get(handles.image_data_options_parent_menu,'children');
                     get(handles.set_resolution_parent_menu_option,'children');
                     handles.hyperparm_quality_report_options;
                     get(handles.hyperparm_quality_report_options,'children')
                     ];
    
    % MenuSelectedFcn is not an available callback in older matlabs, so accomodate that
    menucallbackname = 'MenuSelectedFcn';
	test_menu_handle = get(generic_menus(1));
    if ~isfield(test_menu_handle,'MenuSelectedFcn')
        menucallbackname = 'Callback';
    end
    
    set(generic_menus,menucallbackname,@(hObject,eventdata)svrlsmgui('UpdateCurrentAnalysis',guidata(hObject),hObject))

    % disable things not ready for prime time...
    set([handles.parent_cache_menu handles.beta_options_menu],'enable','off' )
    
    %% Populate the progress bar.
    curvature = .2;
    handles.progressaxes_rectangle_background = rectangle(handles.progressaxes,'Position',[0 0 100 1],'Curvature',curvature,'facecolor',[.95 .95 .9]);
    initial_prog=0;
    handles.progressaxes_rectangle = rectangle(handles.progressaxes,'Position',[0 0 initial_prog 1],'Curvature',curvature,'facecolor',[.8 .9 .7],'tag','progress_rectangle');
    set(handles.progressaxes,'xlim',[0 100],'ylim',[0 1],'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[],'box','off');
    initial_text=''; % nothing.
    handles.progressaxes_text = text(diff(get(handles.progressaxes,'xlim'))/2,diff(get(handles.progressaxes,'ylim'))/2,initial_text, ...
        'horizontalalignment','center','verticalalignment','middle','parent',handles.progressaxes,'tag','progress_text');
    
    %% Try to put image on the interrupt button...
    functions_dir = fileparts(which(mfilename));
    try
        stopicon = fullfile(functions_dir,'other','stop.png');
        x=imresize(imread(stopicon), [22 22]);
        set(handles.interrupt_button,'CData',x)
    catch
        warndlg('Incomplete installation? Not able to locate some of the helper files for SVRLSMGUI. Some aspects of this software may not function.')
    end
    
    %% Try to put image on the open buttons...
    try
        openicon = fullfile(functions_dir,'other','external-link-symbol2.png');
        x=repmat(imresize(imread(openicon), [16 16]),[1 1 3]);
        set([handles.open_output_folder_button handles.open_lesion_folder_button handles.open_score_file_button],'CData',x)
        % todo: make the white backgrounds transparent...
    end