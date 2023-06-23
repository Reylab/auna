function varargout = auna(varargin)
% auna MATLAB code for auna.fig
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @auna_OpeningFcn, ...
                   'gui_OutputFcn',  @auna_OutputFcn, ...
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


% --- Executes just before auna is made visible.
function auna_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% Choose default command line output for auna
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes auna wait for user response (see UIRESUME)
% uiwait(handles.auna);

function auna_CloseRequestFcn(hObject, eventdata, handles)
stop_av(handles);
delete(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = auna_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
hFig = findobj('Tag', 'auna');
set(hFig, 'units','normalized','outerposition',[0 0 1 1]);
set(hFig, 'Visible','off');
imshow('reylab_icon_auna.png','Parent',handles.icon_ax);
disableDefaultInteractivity(handles.icon_ax)
disp('Loading...')
drawnow
start_av(hObject);
drawnow
disp('Done')
set(hFig, 'Visible','on');

% --- Executes on button press in pause_tb.
function pause_tb_Callback(hObject, eventdata, handles)
set(hObject,'Enable','off');
button_state = get(hObject,'Value');
if button_state == 1
    pause(handles.player)
    stop(handles.timer_plot_loop);
else
    resume(handles.player)
    set(handles.add_mark_pb,'Enable','off');
    set(handles.add_start_pb,'Enable','off');
    set(handles.add_end_pb,'Enable','off');

    if  strcmp(get(handles.timer_plot_loop,'Running'),'off')
        set(hObject,'Enable','on');
        start(handles.timer_plot_loop);
    end
end
set(hObject,'Enable','on');


function stop_av(handles)
    stop(handles.timer_plot_loop);
    stop(handles.player)

    
function plot_loop(obj, event,hObject)
    handles = guidata(hObject);
    sr = handles.sr;
    FRlength = handles.par.frame_len;
    new_chunck = (handles.player.CurrentSample + handles.audio_ref) >= sr*FRlength*handles.current_frame ;

    if  new_chunck
        setappdata(hObject,'keyPres',1);
        handles.selected_events = struct();
        cla(handles.multimedia_plot);
        handles = replot_zones(handles);
        handles.current_frame = handles.current_frame + 1;
        ind_beg = (floor(handles.current_frame-1)*FRlength*sr+1);
        ind_end = ceil(handles.current_frame*FRlength*sr);
        if ind_end > handles.lts
             ind_end = handles.lts;
        end
        
        plot(handles.multimedia_plot,handles.ejex_temp(ind_beg:3:ind_end),handles.audio(ind_beg:3:ind_end)); %isn't necesary too much precision
        samp_2_draw = handles.player.CurrentSample+ handles.audio_ref;
        curr_time = handles.ejex_temp(samp_2_draw);

        hold(handles.multimedia_plot,'on')
        borders = xlim(handles.multimedia_plot);
        gap = borders(2)-borders(1);
        if gap >FRlength
            gap = FRlength;
        end
        xlim(handles.multimedia_plot,[handles.ejex_temp(ind_beg) handles.ejex_temp(ind_beg)+gap]);
        ylim(handles.multimedia_plot,[-1 1])
        
        xlabel(handles.multimedia_plot,'Time [s]');
        if ~isempty(handles.concepts)
            concepts_names = fieldnames(handles.concepts);
            for cni=1:numel(concepts_names)
                events = handles.concepts.(concepts_names{cni});
                if ~isempty(events)
                    events = events((events<=handles.ejex_temp(ind_end)*1e3) & (events>=handles.ejex_temp(ind_beg)*1e3));
                    for ie = 1:length(events)
                         e = events(ie);
                         plot(handles.multimedia_plot,[1  1]*e/1e3,[-1 1],'-.','linewidth',2,'color',handles.colors(mod(cni,size(handles.colors,1)),:));
                         text(e/1e3,1,[num2str(cni) ' '],'Parent',handles.multimedia_plot,'Color',handles.colors(mod(cni,size(handles.colors,1)),:),'FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
                    end 
                end
            end    
        end
        handles.lineav = plot(handles.multimedia_plot,[curr_time curr_time],[-1 1],'linewidth',1,'color','k');
        if handles.par.show_fr || handles.par.show_raster
            hold(handles.fr_axes,'off')
            handles.linefr = plot(handles.fr_axes,[curr_time curr_time],[0 handles.plot_counter],'linewidth',1,'color','k');
            hold(handles.fr_axes,'on')
            
            if ~isempty(handles.concepts)
                concepts_names = fieldnames(handles.concepts);
                for cni=1:numel(concepts_names)
                    events = handles.concepts.(concepts_names{cni});
                    if ~isempty(events)
                        events = events((events<=handles.ejex_temp(ind_end)*1e3) & (events>=handles.ejex_temp(ind_beg)*1e3));
                        for ie = 1:length(events)
                             e = events(ie);
                             plot(handles.fr_axes,[1  1]*e/1e3,[0 handles.plot_counter],':','linewidth',2,'color',handles.colors(mod(cni,size(handles.colors,1)),:));
                             text(e/1e3,handles.plot_counter,[num2str(cni) ' '],'Parent',handles.fr_axes,'Color',handles.colors(mod(cni,size(handles.colors,1)),:),'FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
                        end 
                    end
                end    
            end
            
            if handles.plot_counter>1
                yborders = 1:handles.plot_counter-1;
                plot(handles.fr_axes,[ones(1, length(yborders))*handles.ejex_temp(ind_beg); ones(1, length(yborders))*handles.ejex_temp(ind_end)],[yborders;yborders],'--','color',[0,0,0]+0.6);
            end
            xlim(handles.fr_axes,[handles.ejex_temp(ind_beg) handles.ejex_temp(ind_beg)+gap]);
            %xlim(handles.fr_axes,[handles.ejex_temp(ind_beg) handles.ejex_temp(ind_end)]);
            set(handles.fr_axes,'YTick',[])
            set(handles.fr_axes,'XTick',[])
            set(handles.fr_axes,'YMinorGrid','on')
            plot_counter = 1;
            for ch_n = 1: length(handles.par.channels)
                classes = handles.par.classes{ch_n};
                if strcmp(classes,'mu')
                    if handles.par.show_fr
                        plot(handles.fr_axes,handles.ejex_temp(ind_beg:2:ind_end),handles.ch_data{ch_n}.fr(ceil(ind_beg/2):ceil(ind_end/2)),'color',handles.colors(mod(plot_counter,size(handles.colors,1)),:),'linewidth',1);
                    end         
                    
                    if handles.par.show_raster
                        sp = handles.ejex_temp(find(handles.ch_data{ch_n}.sp_index(ind_beg:ind_end))+ind_beg);
                        plot(handles.fr_axes,[sp; sp], [ones(1, length(sp))*(plot_counter-1); ones(1, length(sp))*(plot_counter-1+0.2)],'color','k','linewidth',2);
                    end
                    plot_counter = plot_counter +1;
                else
                    if iscell(classes)
                        classes = cell2mat(classes); %fix when indexing handles.par.classes output cell
                    end
                    for cl_n = 1:length(classes)
                        if handles.par.show_fr
                            plot(handles.fr_axes,handles.ejex_temp(ind_beg:2:ind_end),handles.ch_data{ch_n}.fr(cl_n,ceil(ind_beg/2):ceil(ind_end/2)),'color',handles.colors(mod(plot_counter,size(handles.colors,1)),:),'linewidth',1);
                        end                    
                        if handles.par.show_raster
                            sp = handles.ejex_temp(find(handles.ch_data{ch_n}.index(cl_n,ind_beg:ind_end))+ind_beg);
                            plot(handles.fr_axes,[sp; sp], [ones(1, length(sp))*(plot_counter-1); ones(1, length(sp))*(plot_counter-1+0.2)],'color','k','linewidth',2);
                        end
                        plot_counter = plot_counter +1;
                    end
                end
            end
            
            if handles.plot_counter ~= 0
                ylim(handles.fr_axes,[0 handles.plot_counter])
            end
        end
        guidata(hObject,handles);
        setappdata(hObject,'keyPres',0);
    else
        newcurrent_time = handles.ejex_temp(handles.player.CurrentSample+handles.audio_ref);
        set(handles.lineav, 'Xdata', [newcurrent_time newcurrent_time])
        borders = xlim(handles.multimedia_plot);
        gap = borders(2)-borders(1);
        if newcurrent_time > borders(2)
            xlim(handles.multimedia_plot,[newcurrent_time newcurrent_time+gap]);
        end
        if handles.par.show_fr || handles.par.show_raster
            if newcurrent_time > borders(2)
                xlim(handles.fr_axes,[newcurrent_time newcurrent_time+gap]);
            end
            set(handles.linefr, 'Xdata', [newcurrent_time newcurrent_time])
        end
    end
    
    if handles.player.CurrentSample ~= 1
        set(handles.curr_time_label,'String',[num2str(handles.par.tbeg+(handles.player.CurrentSample + handles.audio_ref - 1)/handles.sr,'%3.3f') 's']);
    elseif handles.current_frame == handles.frame_max
        set(handles.curr_time_label,'String','End');
        stop(handles.player)
    end

    if handles.par.show_mav
        
        current_sample = ceil(handles.player.CurrentSample+handles.audio_ref);

        
        cc = 1;
        for ch_n = 1: length(handles.par.channels)
            ch = handles.par.channels(ch_n);
            classes = handles.par.classes{ch_n};
            if strcmp(classes,'mu')
                aux = handles.ch_data{ch_n}.mav(current_sample);
                handles.clabels{cc,2}.String = sprintf('%s w:%0.1f',handles.clabels{cc,1},aux);                cc = cc +1;
            else
                if iscell(classes)
                    classes = cell2mat(classes); %fix when indexing handles.par.classes output cell
                end
                for cl_n = 1:length(classes)
                    aux = handles.ch_data{ch_n}.mav(cl_n,current_sample);
                    handles.clabels{cc,2}.String = sprintf('%s w:%0.1f',handles.clabels{cc,1},aux);
                    cc = cc +1;
                end
            end
        end
    end
drawnow    

    
    
        
function start_av(hObject)
    clear functions % reset functions, force to reload set_parameters next
    handles = guidata(hObject);
    set(handles.pause_tb,'Value',1)
    handles.par = par_auna();
    handles.current_frame = 0;
    max_channels = length(handles.par.channels);
    plot_counter = 0;
    handles.selected_events = struct();
    setappdata(hObject,'keyPres',0);
    xlim(handles.multimedia_plot,[0 handles.par.frame_len]);
    xlim(handles.multimedia_plot,'manual')
    ylim(handles.multimedia_plot,[-1 1])
           
    
    if isempty(handles.par.folder_base)
        handles.par.folder_base = pwd;
    end
    inds_sep = strfind(handles.par.folder_base,filesep);
    handles.session = handles.par.folder_base(inds_sep(end)+1:end);
    colors = csvread([fileparts(mfilename('fullpath')) filesep 'auna_colors.csv']);    
    handles.colors = colors(2:end,:);
    set(handles.add_mark_pb,'Enable','off');
    set(handles.add_start_pb,'Enable','off');
    set(handles.add_end_pb,'Enable','off');
    if exist(handles.par.data_file,'file')
        load(handles.par.data_file, 'concepts','zones')
        handles.concepts = concepts;
        handles.zones = zones;
    else
        handles.concepts = struct();
        handles.zones = struct();
        disp('Events data not found, start without it');
    end
    
    zone_names = fieldnames(handles.zones);
    concepts_names = fieldnames(handles.concepts);
    listboxItems = cell(length(concepts_names),1);
    for k = 1 : length(concepts_names) 
        listboxItems{k} = ['(' num2str(k) ') ' concepts_names{k} ];
    end
    if isempty(listboxItems)
        set(handles.concepts_pm, 'String', ' ');
        set(handles.concepts_pm,'Enable','off');
    else
        set(handles.concepts_pm,'Enable','on');
        set(handles.concepts_pm, 'String', listboxItems);
    end
    
    listboxItems_zones_pm = cell(length(zone_names),1);
    for k = 1 : length(zone_names) 
        listboxItems_zones_pm{k} = zone_names{k};
    end
    if isempty(listboxItems_zones_pm)
        set(handles.zones_pm, 'String', ' ');
        set(handles.zones_pm,'Enable','off');
    else
        set(handles.zones_pm,'Enable','on');
        set(handles.zones_pm, 'String', listboxItems_zones_pm);
    end
    
    handles.zone_gelements = struct;
    for zi = 1:length(zone_names)
        handles.zone_gelements.(zone_names{zi}) = {};
    end
    handles = replot_zones(handles);
    ch_outnames = cell(max_channels,1);
    if exist(fullfile(handles.par.folder_base,'NSx.mat'),'file')
        lts = [];
        sr = [];
        if isnumeric(handles.par.audio_channel)
            load(fullfile(handles.par.folder_base,'NSx.mat'),'NSx','freq_priority');
            selected = arrayfun(@(x) (x.chan_ID==handles.par.audio_channel)*(find(freq_priority(end:-1:1)==x.sr)),NSx);
            if sum(selected)==0
                error('channel not found in NSx.mat')
            elseif length(nonzeros(selected))>1
                [posch,~] = max(selected);
            else
                posch = find(selected);
            end
            sr(end+1) = NSx(posch).sr;
            lts(end+1) = NSx(posch).lts;
            audio_channel_file = [NSx(posch).output_name NSx(posch).ext];
        else
            sr(end+1) = handles.par.audio_sr;
            aux_info = whos('-file',handles.par.audio_file,handles.par.audio_variable);
            lts(end+1) = max(aux_info.size);
        end

        for ch_n = 1: max_channels
            if exist('freq_priority','var') == false
                load(fullfile(handles.par.folder_base,'NSx.mat'),'NSx','freq_priority');
            end                
            ch = handles.par.channels(ch_n);
            selected = arrayfun(@(x) (x.chan_ID==ch)*(find(freq_priority(end:-1:1)==x.sr)),NSx);
            if sum(selected)==0
                error('channel not found in NSx.mat')
            elseif length(nonzeros(selected))>1
                [posch,~] = max(selected);
            else
                posch = find(selected);
            end
            sr(end+1) = NSx(posch).sr;
            lts(end+1) = NSx(posch).lts;
            ch_outnames{ch_n}= NSx(posch).output_name;
        end

        lts = min(lts);
        sr = unique(sr);
        if length(sr)>1  
            error('Multiple sampling rates')
        end
    else
        audio_channel_file = sprintf('NSX%d.NC5', handles.par.audio_channel);
        load(fullfile(handles.par.folder_base,'NSX_TimeStamps.mat'),'sr','lts');
        for ch_n = 1: max_channels
            ch_outnames{ch_n}=['NSX' num2str(handles.par.channels(ch_n))];
        end
    end

    handles.sr= sr;
    if (sr * handles.par.tbeg) >= lts
        error('tbeg out of recording')
    end    
    if strcmp(handles.par.rec_length,'end') % rec_length in seconds. to decide from which point I playback the sound and show the FR
        tend = lts/sr;
        handles.par.rec_length = tend - handles.par.tbeg;
    else
        tend = handles.par.rec_length + handles.par.tbeg;
    end
    lts = floor(handles.par.rec_length*sr);

    handles.lts = lts;
    handles.frame_max = ceil(handles.par.rec_length/handles.par.frame_len);
    handles.ejex_temp = linspace(handles.par.tbeg,tend,lts);
    
    if handles.par.tbeg==0
        handles.min_record=1;
    else
        handles.min_record = ceil(sr * handles.par.tbeg);
    end
    
    if handles.par.show_fr || handles.par.show_raster || handles.par.show_mav
        ylim(handles.fr_axes,'manual');
        xlim(handles.fr_axes,'manual');
        handles.ch_data = cell(length(handles.par.channels),1);
        if handles.par.show_fr
            half_width_gauss = handles.par.alpha_gauss * handles.par.sigma_gauss;
            sample_period = 1000/sr; % sample period for the spike list - window convolution in ms/sample
            N_gauss = 2*round(half_width_gauss/sample_period)+1; % Number of points of the gaussian window
            int_window = gausswin(N_gauss, handles.par.alpha_gauss);
            int_window = 1000*int_window/sum(int_window)/sample_period;
        end
        handles.ch_data={};
        for ch_n = 1: max_channels
            ch = handles.par.channels(ch_n);
            if strcmp(handles.par.neuroformat,'WC')
                load(fullfile(handles.par.folder_base,sprintf('times_%s.mat',ch_outnames{ch_n})),'cluster_class');
            elseif strcmp(handles.par.neuroformat,'NWB')
                nwb = nwbRead(sprintf('chan_%s.nwb',ch_outnames{ch_n}));
            else
                error('Unknown file format ')
            end
            if strcmp(handles.par.classes{ch_n},'all')
                if strcmp(handles.par.neuroformat,'NWB')
                    handles.par.classes{ch_n} = unique(cluster_class(:,1))';
                else
                    handles.par.classes{ch_n} = nwb.units.id.data;
                end
            end
            classes = handles.par.classes{ch_n};
            if strcmp(classes,'mu')
                 if strcmp(handles.par.neuroformat,'NWB')
                     sorted_times = [];
                     for cii = 1:length(nwb.units.id.data)
                        sorted_times = [sorted_times ;nwb.units.getRow(cii).spike_times{1}];
                     end
                 else
                    sorted_times = cluster_class(:,2);
                 end
                 
                 sorted_times = sorted_times(sorted_times>(handles.par.tbeg*1000))-handles.par.tbeg*1000;
                 handles.ch_data{ch_n}.sp_index = zeros(lts,1);
                 sorted_samples = ceil(sorted_times*sr/1e3);
                 sorted_samples = sorted_samples(sorted_samples<=lts);
                 handles.ch_data{ch_n}.sp_index(sorted_samples) =1;
                 if handles.par.show_fr
                     n_spike_timeline = length(handles.ch_data{ch_n}.sp_index);
                     integ_timeline_stim = conv(handles.ch_data{ch_n}.sp_index, int_window);
                     handles.ch_data{ch_n}.fr = integ_timeline_stim(round(half_width_gauss/sample_period)+1:2:n_spike_timeline+round(half_width_gauss/sample_period));
                     handles.ch_data{ch_n}.frmax = prctile(handles.ch_data{ch_n}.fr,99);
                     handles.ch_data{ch_n}.fr = single(handles.ch_data{ch_n}.fr/handles.ch_data{ch_n}.frmax+(plot_counter)); % rescale 0 to 1, and add offset for plotting
                 end
                 if handles.par.show_mav
                    mav = conv(handles.ch_data{ch_n}.sp_index, ones(1,ceil(handles.par.window_len_mav*handles.sr)))/handles.par.window_len_mav;
                    handles.ch_data{ch_n}.mav = single(mav); % rescale 0 to 1, and add offset for plotting
                 end
                 handles.ch_data{ch_n}.sp_index = uint8(handles.ch_data{ch_n}.sp_index);
                 plot_counter = plot_counter+1;
            else
                max_cls = length(classes);
                handles.ch_data{ch_n}.index = zeros(max_cls,lts);
                handles.ch_data{ch_n}.fr = zeros(max_cls,ceil(lts/2));  %decimate 2
                handles.ch_data{ch_n}.frmax =zeros(max_cls,1);
                if handles.par.show_mav
                    handles.ch_data{ch_n}.mav = zeros(max_cls,lts);
                end
                for cl_n = 1:max_cls
                    class =  classes(cl_n);
                    if strcmp(handles.par.neuroformat,'NWB')
                        sorted_times = nwb.units.getRow(find(nwb.units.id.data==class)).spike_times{1}*1000; %to ms
                    else
                        sorted_times = cluster_class(cluster_class(:,1)==class,2);
                    end
                    
                    sorted_times = sorted_times(sorted_times>(handles.par.tbeg*1000))-handles.par.tbeg*1000;
                    sorted_samples = ceil(sorted_times*sr/1e3);
                    sorted_samples = sorted_samples(sorted_samples<=lts);
                    handles.ch_data{ch_n}.index(cl_n,ceil(sorted_times*sr/1e3)) = 1;
                    if handles.par.show_fr
                        n_spike_timeline = length(handles.ch_data{ch_n}.index(cl_n,:));
                        integ_timeline_stim = conv(handles.ch_data{ch_n}.index(cl_n,:), int_window);
                        fr = integ_timeline_stim(round(half_width_gauss/sample_period)+1:2:n_spike_timeline+round(half_width_gauss/sample_period));
                        handles.ch_data{ch_n}.frmax(cl_n) = max(fr);
                        handles.ch_data{ch_n}.fr(cl_n,:) = single(fr /handles.ch_data{ch_n}.frmax(cl_n)+(plot_counter)); %decimate 2
                    end
                    if handles.par.show_mav
                        mav = conv(handles.ch_data{ch_n}.index(cl_n,:), ones(1,ceil(handles.par.window_len_mav*handles.sr)))/handles.par.window_len_mav;
                        handles.ch_data{ch_n}.mav(cl_n,:) = single(mav(1:lts)); % rescale 0 to 1, and add offset for plotting
                    end
                    plot_counter = plot_counter+1;
                end
            end
        end
    else
        set(handles.fr_axes,'Visible','off')
    end
    
    handles.plot_counter = plot_counter;
    if plot_counter ~= 0
        
        ylim(handles.fr_axes,[0 plot_counter]);
        ylim(handles.axes_labels,[0 plot_counter]);
        xlim(handles.axes_labels,[-1 1])
        yborders = 0:plot_counter;
        plot(handles.axes_labels,[ones(1, length(yborders))*-1; ones(1, length(yborders))],[yborders;yborders],'-','color',[0,0,0]+0.6);
        text_n = 1;
        handles.clabels =cell(0,0);
        for ch_n = 1: max_channels
            ch = handles.par.channels(ch_n);
            classes = handles.par.classes{ch_n};
            
            if strcmp(classes,'mu')
                label = sprintf('Ch: %d \nmu',handles.par.channels(ch_n));
                handles.clabels{end+1,1} = label;
                handles.clabels{end,2} = text(0.1,-0.5+text_n,label,'Parent',handles.axes_labels, 'FontSize',11,'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
                if handles.par.show_fr
                    label = num2str(handles.ch_data{ch_n}.frmax,'%3.2f');
                    text(0.5,text_n-0.08,label,'Parent',handles.axes_labels, 'FontSize',9,'HorizontalAlignment','right','VerticalAlignment','cap')
                end
                text_n = text_n + 1;
            else
                for cl_n = 1:length(classes)
                    label = sprintf('Ch: %d \nCl: %d',handles.par.channels(ch_n),classes(cl_n));
                    handles.clabels{end+1,1} = label;
                    handles.clabels{end,2} = text(0.1,-0.5+text_n,label,'Parent',handles.axes_labels, 'FontSize',11,'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
                    if handles.par.show_fr
                        label = num2str(handles.ch_data{ch_n}.frmax(cl_n),'%3.2f');
                        text(0.5,text_n-0.08,label,'Parent',handles.axes_labels, 'FontSize',8,'HorizontalAlignment','right','VerticalAlignment','cap')
                    end
                    text_n = text_n + 1;
                end
            end
        end
        set(handles.axes_labels,'Visible','off')
    end

    set(handles.ref_label,'String','Current time:');
    handles.ev_time = -1 ;

    handles.audio_ref = 0;

    samples_2_play = floor(sr * handles.par.rec_length);

    if isnumeric(handles.par.audio_channel)
        f1 = fopen(fullfile(handles.par.folder_base,audio_channel_file),'r','l');
        fseek(f1,(handles.min_record-1)*2,'bof');
        y = fread(f1,samples_2_play,'int16=>double');
        fclose(f1);
    else
        y = load(handles.par.audio_file,handles.par.audio_variable);
        y = y.(handles.par.audio_variable);
        y = reshape(y(handles.min_record:end),1,[]);
    end
    ymin = prctile(y,0.5);
    ymax = prctile(y,99.5);
    handles.audio = (y(:)-ymin)/(ymax-ymin)*2-1;
    handles.player = audioplayer(handles.audio(1:end),sr,16);
    timer_period = 0.05;
    play(handles.player);

    handles.timer_plot_loop = timer('Name','plot_loop','TimerFcn',{@plot_loop,hObject},'Period',timer_period,'ExecutionMode','fixedRate');
    guidata(hObject, handles);
    start(handles.timer_plot_loop)
    pause_tb_Callback(handles.pause_tb,[],handles)

% --- Executes on button press in selec_pb.
function selec_pb_Callback(hObject, eventdata, handles)
% hObject    handle to selec_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
uitoolbar = findall(handles.uitoolbar1, '-property', 'Enable');
set( uitoolbar,'Enable', 'off')
if get(handles.pause_tb,'Value') == 0
    set(handles.pause_tb,'Value',1.0);  
    pause_tb_Callback(handles.pause_tb,[],handles);
end
axes(handles.multimedia_plot)
[time_selected, aux,button] = ginput(1);                                          %gets the mouse input
if button == 3
	set(hObject,'Enable','on');
    set( uitoolbar,'Enable', 'on')
    return
end
%set(handles.lineav, 'Xdata', [time_selected time_selected])
sr = handles.sr;
FRlength = handles.par.frame_len;

stop(handles.player)
handles.ev_time = time_selected;
new_init = floor(time_selected*sr);
handles.player = audioplayer(handles.audio(new_init-handles.min_record:end),sr,16);
handles.audio_ref = new_init-handles.min_record;
guidata(hObject, handles);
plot_loop([], [],hObject);
set(handles.add_start_pb,'Enable','on');
set(handles.add_end_pb,'Enable','on');
set(handles.add_mark_pb,'Enable','on');
set(hObject,'Enable','on');    
set( uitoolbar,'Enable', 'on')



function auna_keypressfcn(hObject, eventdata, handles)
 % determine the key that was pressed 
%     setappdata(hObject,'keyPres',0);
 flag = getappdata(hObject,'keyPres');
 if flag == 0
     setappdata(hObject,'keyPres',1); 
     keyPressed = eventdata.Key; 
     switch keyPressed
         case 'subtract'
             uipushtool3_ClickedCallback(hObject, eventdata, handles)
         case 'space'
            button_state = get(handles.pause_tb,'Value');
            if button_state == 1
                set(handles.pause_tb,'Value',0.0);
            else
                set(handles.pause_tb,'Value',1.0);    
            end
            pause_tb_Callback(handles.pause_tb,[],handles);
         case 'leftarrow'
             prev_pb_Callback(handles.prev_pb,[],handles);
         case 'rightarrow'
            next_pb_Callback(handles.next_pb,[],handles);
         case 'return'
             if strcmp(get(handles.add_mark_pb,'Enable'),'on')
                add_mark_pb_Callback(hObject, eventdata, handles)
             end
         otherwise
              keyPressed = keyPressed(end); %for read crrectly the numpad numbers
              if isstrprop(keyPressed,'digit')
                nevent = str2num(keyPressed);
                if nevent>0 && nevent <= length(get(handles.concepts_pm,'String'))
                    set(handles.concepts_pm,'Value',nevent);
                end
              end
           
     end
 end


% --- Executes on button press in save_pb.
function save_pb_Callback(hObject, eventdata, handles)
colors = handles.colors;
concepts = handles.concepts;
zones = handles.zones;
save(handles.par.data_file,'concepts','zones','colors')

function set_param_pb_Callback(hObject, eventdata, handles)
% hObject    handle to set_param_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in set_param_button.
restart = set_par_ui();
if restart
    stop_av(handles);
    start_av(hObject);
end

% --- Executes on button press in prev_pb.
function prev_pb_Callback(hObject, eventdata, handles)
stop(handles.timer_plot_loop);
sr = handles.sr;
FRlength = handles.par.frame_len;
button_state = get(handles.pause_tb,'Value');


stop(handles.player)
if handles.current_frame == 1
    new_init = sr*FRlength*(handles.current_frame-1)+1;
    handles.current_frame = handles.current_frame-1;
    handles.player = audioplayer(handles.audio(new_init:end),sr,16);
    handles.audio_ref = new_init;
    if button_state == 1
        play(handles.player);
        pause(handles.player);
        guidata(hObject, handles);
        plot_loop([],[],hObject);
    else
        play(handles.player);
        guidata(hObject, handles);
        start(handles.timer_plot_loop);
    end
    return
end
new_init = sr*FRlength*(handles.current_frame-2)+1;
handles.current_frame = handles.current_frame-2;
handles.player = audioplayer(handles.audio(new_init:end),sr,16);
handles.audio_ref = new_init;
if button_state == 1
    play(handles.player)
    pause(handles.player);
    guidata(hObject, handles);
    plot_loop([],[],hObject);
else
    play(handles.player);
    guidata(hObject, handles);
    start(handles.timer_plot_loop);
end

set(handles.add_start_pb,'Enable','off');
set(handles.add_end_pb,'Enable','off');
set(handles.add_mark_pb,'Enable','off');

% --- Executes on button press in next_pb.
function next_pb_Callback(hObject, eventdata, handles)
if handles.current_frame >= handles.frame_max
    return
end
stop(handles.timer_plot_loop);
sr = handles.sr;
FRlength = handles.par.frame_len;
button_state = get(handles.pause_tb,'Value');
stop(handles.player);
new_init = sr*FRlength*handles.current_frame+1;
handles.player = audioplayer(handles.audio(new_init:end),sr,16);
handles.audio_ref = new_init;
if button_state == 1
    play(handles.player)
    pause(handles.player)
    guidata(hObject, handles);
    plot_loop([], [],hObject)
    pause_tb_Callback(handles.pause_tb,[],handles);
else
    play(handles.player);
    guidata(hObject, handles);
    start(handles.timer_plot_loop);
end
set(handles.add_start_pb,'Enable','off');
set(handles.add_end_pb,'Enable','off');
set(handles.add_mark_pb,'Enable','off');


function add_mark_pb_Callback(hObject, eventdata, handles)

nevent = get(handles.concepts_pm,'Value');
concepts_names = fieldnames(handles.concepts);
concept_name = concepts_names{nevent};
handles.concepts.(concept_name)(end+1) = handles.ev_time*1e3; % in ms from the beginning of the recording
plot(handles.multimedia_plot,[1  1]*handles.ev_time,[-1 1],'-.','linewidth',2,'color',handles.colors(mod(nevent,size(handles.colors,1)),:));
text(handles.ev_time,1,[num2str(nevent) ' '],'Parent',handles.multimedia_plot,'Color',handles.colors(mod(nevent,size(handles.colors,1)),:),'FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
set(handles.add_mark_pb,'Enable','off');
set(handles.add_start_pb,'Enable','off');
set(handles.add_end_pb,'Enable','off');
guidata(hObject, handles);


% --- Executes on button press in select_event_pb.
function select_event_pb_Callback(hObject, eventdata, handles)
% hObject    handle to select_event_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uitoolbar = findall(handles.uitoolbar1, '-property', 'Enable');
set( uitoolbar,'Enable', 'off')
set(handles.pause_tb,'Value',1.0);    
pause_tb_Callback(handles.pause_tb,[],handles);
rect = getrect(handles.multimedia_plot);
if rect(3)==0 || rect(4) == 0
    if exist('handles.rectan','var')
        delete(handles.rectan);
    end
    handles.selected_events = struct();
    set(uitoolbar,'Enable','on')
    return
end
handles.rectan = plot(handles.multimedia_plot,[rect(1) ,rect(3)+rect(1) ,rect(3)+rect(1),rect(1),rect(1)],[ rect(2),rect(2), rect(4)+rect(2), rect(4)+rect(2),rect(2)],':r','linewidth',2);

etind = rect(1)*1e3;
etend = (rect(1) + rect(3))*1e3;

if ~isempty(handles.concept)
    handles.selected_events = struct();
    concepts_names = fieldnames(handles.concepts);
    for cni=1:numel(concepts_names)
        events = handles.concepts.(concepts_names(cni));
        events = (events>= etind) & (events<= etend);
        if any(events)
            handles.selected_events.(concepts_names(cni))=events;
        end
    end    
end
guidata(hObject, handles);
set(uitoolbar,'Enable','on')

% --- Executes on button press in delete_pb.
function delete_pb_Callback(hObject, eventdata, handles)
% hObject    handle to delete_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.selected_events)
    pause_value = get(handles.pause_tb,'Value');
    if pause_value == 0
        set(handles.pause_tb,'Value',1.0);    
        pause_tb_Callback(handles.pause_tb,[],handles);
    end
    delete(handles.rectan);
    drawnow
    
    concets2rm = fieldnames(handles.selected_events);
    
    %plot white on top of current plots (find and remove them could be
    %easyer) then remove them
    allconcepts = fieldnames(handles.concepts);
    for ic=1:length(concets2rm)
        deev = find(handles.selected_events.(concets2rm(ic)));
        nevent = find(cellfun(@(x) strcmp(x,concets2rm(ic)),allconcepts));
        for i = 1:length(deev)
            tevent = handles.concets.(concets2rm{ic})(deev(i))/1e3;
            plot(handles.multimedia_plot,[1  1]*tevent,[-1 1],'-.','linewidth',2,'color','w');
            text(tevent,1,[num2str(nevent) ' '],'Parent',handles.multimedia_plot,'Color','w','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
            if handles.par.show_fr || handles.par.show_raster
             plot(handles.fr_axes,[1  1]*tevent,[0 handles.plot_counter],':','linewidth',2,'color','w');
             text(tevent,handles.plot_counter,[num2str(nevent) ' '],'Parent',handles.fr_axes,'Color','w','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top');
            end
        end
    	handles.concets.(concets2rm{ic})(handles.selected_events.(concets2rm(ic))) = [];
    end
    
    handles.selected_events = struct();
    guidata(hObject, handles)
    if pause_value == 0
        set(handles.pause_tb,'Value',0);    
        pause_tb_Callback(handles.pause_tb,[],handles);
    end
end
    
% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
sr = handles.sr;
FRlength = handles.par.frame_len;
ind_beg = (floor(handles.current_frame-1)*FRlength*sr+1);
ind_end = ceil(handles.current_frame*FRlength*sr);
if ind_end > handles.lts
	ind_end = handles.lts;
end

xlim(handles.multimedia_plot,[handles.ejex_temp(ind_beg) handles.ejex_temp(ind_end)]);
ylim(handles.multimedia_plot,[-1 1])
if handles.par.show_fr || handles.par.show_raster
    xlim(handles.fr_axes,[handles.ejex_temp(ind_beg) handles.ejex_temp(ind_end)]);
end
% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject, eventdata, handles)

pause_value = get(handles.pause_tb,'Value');

if pause_value == 0
    set(handles.pause_tb,'Value',1.0);    
    pause_tb_Callback(handles.pause_tb,[],handles);
end
rect = getrect(handles.multimedia_plot);
if rect(3)==0 || rect(4) == 0
    return
end

tind = rect(1);
tend = (rect(1) + rect(3));
xlim(handles.multimedia_plot,[tind tend]);
if handles.par.show_fr || handles.par.show_raster
    xlim(handles.fr_axes,[tind tend]);
end

if pause_value == 0
    set(handles.pause_tb,'Value',0.0);    
    pause_tb_Callback(handles.pause_tb,[],handles);
end

function auna_WindowKeyReleaseFcn(hObject, eventdata, handles)
    setappdata(hObject,'keyPres',0); 

% --- Executes during object creation, after setting all properties.
function time_input_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function multimedia_plot_CreateFcn(hObject, eventdata, handles)
	hObject.Toolbar.Visible = 'off';
    
function fr_axes_CreateFcn(hObject, eventdata, handles)
	hObject.Toolbar.Visible = 'off';
    
% --- Executes on button press in go2time.
function time_input_Callback(hObject, eventdata, handles)
    jump = str2num(handles.time_input.String);
    if isempty(jump)
        set(handles.time_input,'String','input time (sec)');
        return
    end

% --- Executes on button press in go2time.
function go2time_Callback(hObject, eventdata, handles)
    jump = str2num(handles.time_input.String);
    if isempty(jump)
        return
    end
    if jump < handles.par.tbeg
        return
    end
    sr = handles.sr;
    FRlength = handles.par.frame_len;
    frame = floor((jump - handles.par.tbeg)/FRlength); %previous frame to update in plot_loop

    if frame >= handles.frame_max
        return
    end

    stop(handles.timer_plot_loop);
    pause_value = get(handles.pause_tb,'Value');
    if pause_value == 0
        set(handles.pause_tb,'Value',1.0);    
        pause_tb_Callback(handles.pause_tb,[],handles);
    end    
    stop(handles.player)

    handles.audio_ref = sr*FRlength*(frame) +1;
    new_init = handles.audio_ref;
    handles.current_frame = frame;
    handles.player = audioplayer(handles.audio(new_init:end),sr,16);
    play(handles.player)
    pause(handles.player);
    guidata(hObject, handles);
    plot_loop([],[],hObject);

    if pause_value == 0
        set(handles.pause_tb,'Value',0);    
        pause_tb_Callback(handles.pause_tb,[],handles);
    end 
    set(handles.add_mark_pb,'Enable','off');
    set(handles.add_start_pb,'Enable','off');
    set(handles.add_end_pb,'Enable','off');

% --- Executes on selection change in zones_pm.
function zones_pm_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function zones_pm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_annotations_pb.
function set_annotations_pb_Callback(hObject, eventdata, handles)
    app = configuration_manager(hObject);
    waitfor(app)
    handles = guidata(hObject);
    handles = replot_zones(handles);
    guidata(hObject, handles);
    
    
% --- Executes on button press in add_start_pb.
function add_start_pb_Callback(hObject, eventdata, handles)
nzone = get(handles.zones_pm,'Value');
zones_names = fieldnames(handles.zones);
zones_name = zones_names{nzone};

handles.zones.(zones_name).limits(1)= handles.ev_time*1e3; % in ms from the beginning of the recording
handles = replot_zones(handles);

set(handles.add_mark_pb,'Enable','off');
set(handles.add_start_pb,'Enable','off');
set(handles.add_end_pb,'Enable','off');
guidata(hObject, handles);


% --- Executes on button press in add_end_pb.
function add_end_pb_Callback(hObject, eventdata, handles)
nzone = get(handles.zones_pm,'Value');
zones_names = fieldnames(handles.zones);
zone_name = zones_names{nzone};

handles.zones.(zone_name).limits(2)= handles.ev_time*1e3; % in ms from the beginning of the recording
handles = replot_zones(handles);

set(handles.add_mark_pb,'Enable','off');
set(handles.add_start_pb,'Enable','off');
set(handles.add_end_pb,'Enable','off');
guidata(hObject, handles);


% --- Executes on button press in create_rasters_pb.
function create_rasters_pb_Callback(hObject, eventdata, handles)
    app = raster_menu(handles);
    waitfor(app)

       
function handles = replot_zones(handles)
gelements_fields = fieldnames(handles.zone_gelements);
for zi = 1:length(gelements_fields)
	cellfun(@(x) delete(x),handles.zone_gelements.(gelements_fields{zi}))
end
handles.zone_gelements = struct;
zones_names = fieldnames(handles.zones);
    
for zi = 1:length(zones_names)
    zone_name = zones_names{zi};
    elements = {};
    if all(~isnan(handles.zones.(zone_name).limits))
        dt = diff(handles.zones.(zone_name).limits)/1e3;
        if dt<0
            error('Zone Start after End')
        end
        r=rectangle(handles.multimedia_plot,'Position',[handles.zones.(zone_name).limits(1)/1e3 -1 dt 2],'FaceColor',[handles.colors(mod(zi,size(handles.colors,1)),:), 0.08], 'LineStyle','none');
        tt1=text(handles.zones.(zone_name).limits(1)/1e3,1,zone_name,'Parent',handles.multimedia_plot,'Color',handles.colors(mod(zi,size(handles.colors,1)),:),'FontSize',8,'HorizontalAlignment','left','VerticalAlignment','top');
        tt2=text(handles.zones.(zone_name).limits(2)/1e3,1,zone_name,'Parent',handles.multimedia_plot,'Color',handles.colors(mod(zi,size(handles.colors,1)),:),'FontSize',8,'HorizontalAlignment','right','VerticalAlignment','top');
        elements = {r, tt1,tt2};
    end
    handles.zone_gelements.(zone_name) = elements;
end
