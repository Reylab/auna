function par=par_av_gui() 
 
par.folder_base = []; % [] for pwd 
par.tbeg = 0;
par.rec_length = 'end'; % or number of seconds from tbeg 
par.frame_len = 20; % max chunk of data ploted in seconds 
par.neuroformat = 'WC'; % WC or NWB 
par.channels = [2076]; %for multichannel select 
par.audio_channel = 2129; %if int nc5 otherwise filename of a .mat to load
par.data_file = 'auna_data.mat';
%parameters for mat audio file
par.audio_file = 'test.mat';
par.audio_sr = 30000; %in Hz
par.audio_variable = 'data';
 
%Firing rate parameters 
par.show_fr = true; 
par.sigma_gauss = 49.42; 
par.alpha_gauss = 3.035; %last value of gaussian 0.01 0.025 
% lenght of gaussian window = 2* alpha_gauss * sigma_gauss 

%Moving average parameters 
par.show_mav = true; 
par.window_len_mav = 20; %in seconds


par.show_raster = true; 


par.classes = {[1,2]}; % or 'mu' for multi unit, 'all' for all classes or cell
                    % of the same length as channels, with the classes for
                    % each.
 
par.nrows = 5;  %number of rows in the matrix of frames in video mode 
par.fr_prev = 25;

%zone parameters
par.default_pre_time_ms = 1000;
par.default_post_time_ms = 1000;
end 
 

 
 
