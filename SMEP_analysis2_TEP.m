%% SMEP - SCRIPT 2: TEP preprocessing & analysis
% written by Dominika for the SMEP project (2023)
% 
% The SMEP project explores the recording of spinal motor evoked potentials
% (SMEP) and their use to chanracterize functional changes along the motor
% pathway. To this end, primary motor cortex was stimulated using
% supra-threshold TMS (MagVenture; 120%rMT) and the evoked activity was
% registered at different levels of the cortico-motor pathway:
%   1) TMS-evoked potentials (TEPs) were recorded using 32-channel EEG system
%   (NeurOne, Bittium + EasyCap mounted with TMS compatible multitrodes)
%   2) spinal activity was recorded with a custom 13-channel spinal grid 
%   (multitrodes) + 4 side channels placed on trapezius muscle
%   3) motor-evoked potentials (MEPs) were recorded form the FDI muscle
% 
% Following script provides the code for complete processing of TEPs
% - each section of the script performs one processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - output variables are gathered in an output structure 'SMEP.mat'
% - when the manual input is needed (ICA, visual inspection), the process
%   is performed in the GUI, then the corresponding section of the
%   script takes care of the logfile update and encodes the output
%   parameters to the ouput structure
% 
% 1) PREPROCESSING
%       - assign electrode coordinates
%       - re-reference to common average
%       - epoch around the event [-1 2]s
%       - remove the DC shift and linearly detrend
%       - interpolate the TMS artifact [-0.006 0.003]s
%       - downsample to 2kHz (1/10)
%       - baseine correct by mean subtraction [-0.25 -0.006]s
%       - remove extra event categories
%       - adds global tags to all trials
% 
% 2) REMOVE MISSED OR FAULTY TRIALS
%       - requires manual input of missed trials: a cell 'missed' 
%           - first column indicates block number
%           - second column contains a vector of indexes of missed trials
% 
% 3) PARSE PER CONDITION
%       - loads the stimulation order (needs to be in the data folder!)
%       - labels the epochs with appropriate eventcodes
%       - merges epochs into three condition datasets 
% 
% 4) SUBSET ACCORDING TO MEPS
%       - creates two new datasets:
%           1) without epochs that were discarded based on increased motor
%           activity at the baseline --> 'motor_bl'
%           2) without epochs that were discarded based on the size of
%           peak-to-peak MEP amplitude --> 'motor_zero' 
%               - M1 TEPs: if the amplitude in the window of interest
%               [0.015, 0.050]s is under 30 microV (threshold +- 15)
%               - CTRL TEPs: if the amplitude any time post stimulus
%               exceeds 30 microV (threshold +- 15)
% 
% 5) EXPORT FOR EEGLAB
%       - saves as a .set file
% 
% 6) EEGLAB: BAD CHANNELS & EPOCHS
%       - starts a new EEGLAB session
%       - loads all datasets
%       - interpolates bad channels, if necessary
%       - lets user manualy identify artifactual epochs and discards them
%       - saves the new dataset --> 'visual'
% 
% 7) EEGLAB: SSP-SIR
%       - calculates SSP filter for individual datasets
%       - performs SSP-SIR on each dataset using all three filters
% 
% 8) EEGLAB: FREQUENCY FILTERS + ICA 
%       - 
%       - uses ICA to remove artifact not timelocked to the stimulus 
%           - ocular artifacts 
%           - TMS-independent muscular activity
%           - electrode noise


%% PARAMETERS
clear all; clc;

% subject
subject = 8;
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% dataset
measure = 'TEP';
block = [1:15];
condition = {'M1_single', 'M1_paired', 'CTRL'}; 

% choose relevant directories
folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave + eeglab masterfiles
folder_data = uigetdir(pwd, 'Choose the data folder');              % processed data
folder_output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file
cd(folder_data)

% load output structure 
load([folder_output '\SMEP.mat']);

% create logfile filename, launch the TEP section
filename = sprintf('%s\\Logfiles\\SMEP S%s %s.txt', folder_output, subj, SMEP.info(subject).date); 
% logfile_entry('heading', filename);

% load the finish sound
% load handel
load gong
soundwave = y; clear y Fs

%% 1) PREPROCESSING
% ----- section input -----
suffix = {'reref' 'ep' 'dc' 'interp' 'ds' 'bl'};
eventcode = 'Stimulation';
param.epoch = [-1 2];
param.interp = [-0.006 0.003];
param.ds_ratio = 10;
param.baseline = [-0.25 -0.006];
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% preprocess in letswave
for b = block
    fprintf('block %d:\n', b)

    % load data and header
    [header, data] = CLW_load(sprintf('S%s %s b%d', subj, measure, b));   

    % rereference to common average
    fprintf('rereferencing to average...')
    [header, data, ~] = RLW_rereference(header, data, 'apply_list', {header.chanlocs(1:30).labels}, 'reference_list', {header.chanlocs(1:30).labels});

    % segment 
    fprintf('epoching...')
    [header, data, ~] = RLW_segmentation(header, data, {{eventcode}}, 'x_start', param.epoch(1), 'x_duration', param.epoch(2) - param.epoch(1));
    
    % remove DC and apply linear detrend
    fprintf('removing DC + detrending...')
    [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);

    % interpolate TMS artifact
    fprintf('interpolating TMS artifact...')
    [header, data, ~] = RLW_suppress_artifact_event(header, data,...
        'xstart', param.interp(1), 'xend',  param.interp(2), 'event_code', eventcode, 'interp_method', 'pchip');

    % downsample
    fprintf('downsampling...')
    [header, data, ~] = RLW_downsample(header, data, 'x_downsample_ratio', param.ds_ratio);

    % baseline correct
    fprintf('subtracting baseline...')
    [header, data, ~] = RLW_baseline(header, data, 'operation', 'subtract', 'xstart', param.baseline(1), 'xend', param.baseline(2));
    
    % add the global event tag
    for h = 1:length(header.events)
        header.global_tags(h, :) = [b, h];
    end

    % save for letswave 
    fprintf('done.\n')
    for s = 1:length(suffix)
        header.name = [suffix{s} ' ' header.name];
    end
    CLW_save([], header, data);
end

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% create prefix
for s = 1:length(suffix)
    if s == 1
        prefix = suffix{s};
    else
        prefix = [suffix{s} ' ' prefix];
    end
end

% assign electrode coordinates
fprintf('Assigning electrode coordinates... ')
for b = block   
    % load the data
    option = struct('filename', sprintf('%s S%s %s b%d.lw6', prefix, subj, measure, b));
    lwdata = FLW_load.get_lwdata(option);
    
    % assign electrode coordinates
    option = struct('filepath', 'F:\letswave7-master\letswave7-master\res\electrodes\spherical_locations\Standard-10-20-Cap81.locs', ...
        'suffix', '', 'is_save', 1);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
end
fprintf('done.\n')

% encode to the logfile 
logfile_entry('preprocessing', filename, 'parameters', param);
clear suffix eventcode param b h s data header lwdata option
sound(soundwave)

%% 2) REMOVE MISSED OR FAULTY TRIALS
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% identify missed epochs
missed = SMEP.MEP(subject).missed_epochs;

% discard missing/faulty epochs if necessary
if length(missed{1}) > 0
    fprintf('discarding missed/faulty epochs...\n')
    
    % cycle through blocks containing missing trials
    for b = 1:size(missed, 1)
        % load the dataset
        option = struct('filename', sprintf('%s S%s %s b%d', prefix, subj, measure, missed{b, 1}));
        lwdata = FLW_load.get_lwdata(option);

        % identify epochs to keep
        epochs2keep = {};
        epochs = 1:size(lwdata.data, 1);
        epochs = epochs(~ismember(epochs, missed{b, 2}));
        if isempty(epochs)
            continue
        else
            for e = 1:length(epochs)
                epochs2keep{e} = num2str(epochs(e));
            end
        end
        
        % update global tags
        lwdata.header.global_tags = lwdata.header.global_tags(epochs, :);

        % discard indicated epochs 
        option = struct('type', 'epoch', 'items', {epochs2keep}, 'suffix', '', 'is_save', 1);
        lwdata = FLW_selection.get_lwdata(lwdata, option);    
    end
    
    % encode to the logfile 
    logfile_entry('missed', filename, 'missed_epochs', missed);
    fprintf('done.\n')
else    
    % encode to the logfile 
    logfile_entry('missednot', filename);
    fprintf('no missing trials found.\n')
end

% save the output structure
save([folder_output '\SMEP.mat'], 'SMEP');
clear missed b option lwdata epochs epochs2keep e

%% 3) PARSE PER CONDITION
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% load order of stimulation --> stim_order
load([folder_output '\TMS stimulation protocols\SMEP_' subj '_stim_order.mat'])

% label and cluster by condition
fprintf('parsing per condition...\n')
counter = 1;
global_tags = {[], [], []};
for b = block
    % get the dataset
    option = struct('filename', sprintf('%s S%s %s b%d.lw6', prefix, subj, measure, b));
    lwdata = FLW_load.get_lwdata(option);
    
    % verify the condition
    for c = 1:length(condition)
        if strcmp(condition{c}, stim_order{b})
            % replace event codes
            for e = 1:length(lwdata.header.events)
                lwdata.header.events(e).code = condition{c};
            end
            
            % rename & save into a dataset for merging
            lwdata.header.name = sprintf('S%s %s', subj, measure);
            data2merge(c, counter) = lwdata;
            
            % append to global tag variable
            global_tags{c} = cat(1, global_tags{c}, lwdata.header.global_tags);
        end
    end
    
    % update counter if necessary
    if mod(b, 3) == 0
        counter = counter + 1;
    end    
end

% merge blocks according to conditions
fprintf('merging datasets...\n')
for c = 1:length(condition)
    % subset the dataset
    lwdataset = data2merge(c, :)';
    
%     % in case a whole block needs to be removed
%     if c == 1
%        lwdataset = lwdataset([1:3, 5]);
%     end

    % merge subjects & save 
    option = struct('type', 'epoch', 'suffix', condition{c}, 'is_save', 1);
    lwdata = FLW_merge.get_lwdata(lwdataset, option); 
end
fprintf('done.\n')

% update global tags  
addpath(genpath([folder_toolbox '\letswave6-master']));
for c = 1:length(condition)
    % load header
    [header, data] = CLW_load(sprintf('%s S%s %s', condition{c}, subj, measure)); 
    
    % replace global tags, if the lengths match
    if header.datasize(1) == length(global_tags{c})
        header.global_tags = global_tags{c};
        CLW_save([], header, data);
    else
        error('Error: number of found epochs does not match the length of the global tags file!')
    end
end

% encode to the logfile 
logfile_entry('merged', filename); 
clear b c e option lwdata folder_input stim_order data2merge counter lwdataset header data global_tags

%% 4) SUBSET ACCORDING TO MEPS
% ----- section input -----
suffix = {'motor_bl' 'motor_zero'};
MEP_interval = [0.015, 0.05];
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% cycle through conditions
fprintf('subsetting based on motor activity: ')
for c = 1:length(condition)
    fprintf('...%s ', condition{c})
    
    % load the data
    option = struct('filename', sprintf('%s S%s %s', condition{c}, subj, measure));
    lwdata = FLW_load.get_lwdata(option);
        
    % idetify epochs with baseline motor activity 
    epochs2ditch = SMEP.MEP(subject).baseline_discarded{c};
    epochs2keep = {};
    epochs = 1:size(lwdata.data, 1);
    epochs = epochs(~ismember(epochs, epochs2ditch));
    if isempty(epochs)
        continue
    else
        for e = 1:length(epochs)
            epochs2keep{e} = num2str(epochs(e));
        end
    end
    
    % update global tags
    lwdata.header.global_tags = lwdata.header.global_tags(epochs, :);

    % discard indicated epochs 
    option = struct('type', 'epoch', 'items', {epochs2keep}, 'suffix', suffix{1}, 'is_save', 1);
    lwdata = FLW_selection.get_lwdata(lwdata, option);  
    
    % idetify epochs with baseline zero MEP (in case of M1 stimulation) or
    % post-stimulus motor activity (in case of CTRL stimulation)
    epochs2ditch = SMEP.MEP(subject).zero_discarded{c};
    epochs2keep = {};
    epochs = 1:size(lwdata.data, 1);
    epochs = epochs(~ismember(epochs, epochs2ditch));
    if isempty(epochs)
        continue
    else
        for e = 1:length(epochs)
            epochs2keep{e} = num2str(epochs(e));
        end
    end
    
    % update global tags
    lwdata.header.global_tags = lwdata.header.global_tags(epochs, :);
    
    % discard indicated epochs 
    option = struct('type', 'epoch', 'items', {epochs2keep}, 'suffix', suffix{2}, 'is_save', 1);
    lwdata = FLW_selection.get_lwdata(lwdata, option); 
    
    % verify that TEP global tags match MEP global tags
    load(sprintf('%s %s %s S%s MEP.lw6', suffix{2}, suffix{1}, condition{c}, subj), '-mat')
    for e = 1:length(lwdata.header.global_tags)
        if lwdata.header.global_tags(e, :) ~= header.global_tags(e, :)
            error('Error: Global tags do not match between TEP and MEP files!')
            break
        end
    end 
end
fprintf('...done.\n')

% read out MEP parameters
param.threshold = SMEP.MEP(subject).baseline_threshold;
param.MEP_interval = MEP_interval;

% encode to the logfile 
logfile_entry('motor_filter', filename, 'parameters', param); 
clear suffix MEP_interval c option lwdata header epochs2ditch epochs2keep epochs e param

%% 5) EXPORT FOR EEGLAB
% ----- section input -----
prefix = {'motor_zero motor_bl' 'motor_bl' 'motor_zero motor_bl'};
suffix = 'preprocessed';
ref = 'averef';
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export in .set format
fprintf('exporting for EEGLAB: ')
for c = 1:length(condition)
    fprintf('...%s ', condition{c})
    
    % load the data
    option = struct('filename', sprintf('%s %s S%s %s', prefix{c}, condition{c}, subj, measure));
    lwdata = FLW_load.get_lwdata(option);
    
    % export in .set format
    name = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix);
    export_EEGLAB(lwdata, name, ref);
end
fprintf('...done.\n')

% encode to the logfile 
logfile_entry('eeglab', filename); 
clear prefix suffix ref c lwdata name

%% 6) EEGLAB: BAD CHANNELS & EPOCHS
% ----- section input -----
suffix = {'preprocessed' 'visual'};
% ------------------------- 
% add eeglab and fastica to the of search path
addpath(fullfile(folder_toolbox,'eeglab2022.1'));
addpath(fullfile(folder_toolbox,'FastICA_25'));

% launch eeglab and generate an empty EEGLAB structure
eeglab 

% cycle through conditions
for c = 1:length(condition)
    % load the dataset
    name = sprintf('%s %s S%s %s.set', measure, condition{c}, subj, suffix{1});  
    EEG = pop_loadset('filename', name, 'filepath', folder_data);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw  
    
    % remove ECG channels 
    EEG = pop_select(EEG, 'nochannel', {'ECGleft','ECGright','ECGdown'});
    EEG.ref = 'averef';
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'overwrite', 'on', 'gui', 'off'); 
    eeglab redraw
    
    % visualize the average response to identify possible bad channels
    figure; 
    pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);
    sgtitle(sprintf('S%s - %s: original data', subj, condition{c}))
    
    % remove bad channels if necessary
    answer = inputdlg('Which channels do you want to interpolate?', 'Bad channels', [1 35], {'none'});
    if strcmp(answer, 'none')
        interpolated{c} = '';
    else
        % identify channels to interpolate
        interpolated{c} = '';
    end
    
    % remove bad epochs
    pop_eegplot(EEG, 1, 1, 1);
    waitfor(gcf); 
    EEG = eeg_checkset(EEG);
    
    % extract information about rejected epoch
    char_start = 'EEG = pop_rejepoch( EEG, [';
    char_end = '] ,0)';
    epochs_start = strfind(EEG.history, char_start) + length(char_start);
    epochs_end = strfind(EEG.history, char_end) - 1;
    eval(['discarded.epochs_rejected{c} = [' EEG.history(epochs_start:epochs_end) '];']);
    discarded.epochs_kept(c) = length(EEG.epoch);
    
    % update global tags
    EEG.global_tags(discarded.epochs_rejected{c}, :) = [];
    
    % save dataset
    name = sprintf('%s %s S%s %s.set', measure, condition{c}, subj, suffix{2}); 
    pop_saveset(EEG, 'filename', name, 'filepath', folder_data);
end

% fill in the output structure and save
SMEP.TEP(subject).ID = subject;
SMEP.TEP(subject).interpolated = interpolated;
SMEP.TEP(subject).visual_discarded = discarded.epochs_rejected;
SMEP.TEP(subject).visual_kept = discarded.epochs_kept;
save([folder_output '\SMEP.mat'], 'SMEP');

% encode to the logfile
logfile_entry('subset', filename, 'interpolated', answer, 'discarded', discarded); 
clear suffix c name answer interpolated discarded

%% 7) EEGLAB: SSP-SIR
% ----- section input -----
suffix = 'sspsir';
time_range = [-6, 50];
% ------------------------- 
% filter individual datasets
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
        
    % visual check - before SSP-SIR
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Original data');
    sgtitle(sprintf('S%s - %s: original data', subj, condition{c}))
    
    % SSP-SIR - spherical model 
    name = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix); 
    [EEG] = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', time_range, 'PC', []);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'setname', name, 'overwrite', 'on', 'gui', 'off'); 
    name = sprintf('%s %s S%s %s.set', measure, condition{c}, subj, suffix); 
    EEG.filename = name;
    eeglab redraw
    EEG = eeg_checkset(EEG);
    
    % visual check - after SSP-SIR
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Filtered data');
    sgtitle(sprintf('S%s - %s: filtered data', subj, condition{c}))
    
    % extract information about rejected components
    param(c) = EEG.SSPSIR;
end

% fill in the output structure and save
SMEP.TEP(subject).sspsir_PC = [param(1).PC, param(2).PC, param(3).PC];
SMEP.TEP(subject).sspsir_params = rmfield(param, 'PC');
save([folder_output '\SMEP.mat'], 'SMEP');

% apply ssp-sir filters of other compared datasets
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
            
    % SSP-SIR 
    other_datasets = [1,2,3];
    other_datasets = other_datasets(other_datasets ~= c);
    for a = other_datasets
        [EEG_trials] = SSP_SIR_trials(EEG, SMEP.TEP(subject).sspsir_params(c).L_ave, ...
            SMEP.TEP(subject).sspsir_params(a).topographies, ...
            SMEP.TEP(subject).sspsir_params(a).filter, []);
        EEG.data =  EEG_trials.data;
    end
    eeglab redraw
    
    % baseline correct 
    EEG = pop_rmbase(EEG, param.baseline, []);
    
    % visual check - after all SSP-SIR rounds
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Final data');
    sgtitle(sprintf('S%s - %s: final data', subj, condition{c}))
    
    % save dataset
    name = sprintf('%s %s S%s %s.set', measure, condition{c}, subj, suffix); 
    pop_saveset(EEG, 'filename', name, 'filepath', folder_data);
end

% encode to the logfile
logfile_entry('ssp-sir', filename, 'PC', SMEP.TEP(subject).sspsir_PC); 
clear suffix time_range c name param EEG.trials

%% 8) EEGLAB: ICA 
% ----- section input -----
suffix = 'ica';
param.baseline = [-0.25 -0.006];
param.compression = 20;
% ------------------------- 
% % load data, if necessary
% eeglab
% for c = 1:length(condition)
%     % load the dataset
%     name = sprintf('%s %s S%s sspsir.set', measure, condition{c}, subj);  
%     EEG = pop_loadset('filename', name, 'filepath', folder_data);
%     [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
%         
%     % baseline correct 
%     EEG = pop_rmbase(EEG, param.baseline, []);
%     eeglab redraw  
%     
%     % modify dataset info
%     EEG.subject = subj; 
%     EEG.session = 1;
% end

% run ICA separately on all datasets 


% run ICA separately for each dataset
for c = 2:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw

    % set the compression level
    EEG = pop_tesa_pcacompress(EEG, 'compVal', param.compression, 'plot', 'off');
    
    % visual check before ICA
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Filtered data');
    sgtitle(sprintf('S%s - %s: frequency-filtered data', subj, condition{c}))
    
    % run ICA 
    EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off');
    EEG = pop_tesa_compplot(EEG, 'figSize', 'large', 'plotTimeX', [-0.2 0.5], 'plotFreqX', [1 100],...
        'freqScale', 'log', 'saveWeights', 'off');
    
    % extract ICA parameters
    SMEP.TEP(subject).ICA(c) = rmfield(EEG.icaCompClass.Manual1, 'name');
    save([folder_output '\SMEP.mat'], 'SMEP');
    
    % baseline correct post ICA
    EEG = pop_rmbase(EEG, param.baseline, []);
    
    % rename
    EEG.setname = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix); 
    
    % visual check after ICA
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Data after ICA');
    sgtitle(sprintf('S%s - %s: data after ICA', subj, condition{c}))
    
    % save dataset
    pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', folder_data);
end

clear param suffix 

%% 8) EEGLAB: INTERPOLATE LEFTOVER ARTIFACT + FREQUENCY FILTER
% ----- section input -----
suffix = 'filtered';
param.filter_butter = [1, 80, 4];
param.filter_notch = [50, 2];
param.baseline = [-0.25 -0.006];
% ------------------------- 
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
    
    % butterworth bandpass filter
    name = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix); 
    EEG = pop_tesa_filtbutter(EEG, param.filter_butter(1), param.filter_butter(2), param.filter_butter(3), 'bandpass');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'setname', name, 'overwrite', 'on', 'gui', 'off'); 
    
    % 50Hz notch filter
    EEG = pop_eegfiltnew(EEG, 'locutoff', param.filter_notch(1) - param.filter_notch(2), ...
        'hicutoff', param.filter_notch(1) + param.filter_notch(2), 'revfilt', 1, 'plotfreqz', 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'overwrite', 'on', 'gui', 'off'); 
    eeglab redraw;
end

%% FUNCTIONS
function logfile_entry(entry, filename, varargin)
    switch entry
        case 'heading'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, 'TMS-EVOKED POTENTIALS\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, '1	dataset imported using matlab script ''SMEP_import.m''\r\n');
            fprintf(fileID, '		- BIN data loaded from NeurOne output\r\n');
            fprintf(fileID, '		- 30 scalp channels (no mastoids) + 3 ECG channels\r\n');
            fprintf(fileID, '		- only ''Stimulation'' category kept\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'preprocessing'
            a = find(strcmpi(varargin, 'parameters'));
            param = varargin{a + 1};
            fileID = fopen(filename, 'a');
            fprintf(fileID, 'following processing steps were performed using matlab script ''SMEP_analysis2_TEP.m''\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '2	preprocessing\r\n');
            fprintf(fileID, '		- default electrode coordinates assigned\r\n');
            fprintf(fileID, '		- re-referenced to common average\r\n');
            fprintf(fileID, sprintf('		- segmentation relative to events [%d %d]s\r\n', param.epoch(1), param.epoch(2)));
            fprintf(fileID, '		- DC shift removal + linear detrend\r\n');
            fprintf(fileID, sprintf('		- TMS artifact interpolated at [%.3f %.3f]s\r\n', param.interp(1), param.interp(2)));
            fprintf(fileID, sprintf('		- downsampled to %dkHz (1/%d)\r\n', 20/param.ds_ratio, param.ds_ratio));
            fprintf(fileID, sprintf('		- baseline corrected by mean subtraction [%.3f %.3f]s\r\n', param.baseline(1), param.baseline(2)));
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'missed'
            a = find(strcmpi(varargin, 'missed_epochs'));
            missed_epochs = varargin{a + 1};
            fileID = fopen(filename, 'a');
            fprintf(fileID, '3	missed /faulty epochs were discarded:\r\n');
            for e = 1:size(missed_epochs, 1)
                fprintf(fileID, sprintf('		- block %d --> trial(s) %s\r\n', missed_epochs{e, 1}, num2str(missed_epochs{e, 2})));
            end
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'missednot'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '3	no epochs were marked as missed or faulty\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
        
        case 'merged'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '4	data were labelled and split to three datasets according to the condition\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'motor_filter'
            a = find(strcmpi(varargin, 'parameters'));
            param = varargin{a + 1};
            fileID = fopen(filename, 'a');
            fprintf(fileID, '5	subset datasets were created according to recorded motor activity:\r\n');
            fprintf(fileID, '		1) without epochs that were discarded based on increased motor activity at the baseline\r\n');
            fprintf(fileID, '		2) without epochs that were discarded based on the size of peak-to-peak MEP amplitude\r\n');
            fprintf(fileID, sprintf('		- M1 TEPs: if the amplitude in the window of interest [%.3f %.3f]s\r\n', param.MEP_interval(1), param.MEP_interval(2)));
            fprintf(fileID, sprintf('		is under 30 microV (threshold +- %d)\r\n', param.threshold));
            fprintf(fileID, sprintf('		- CTRL TEPs: if the amplitude any time post stimulus exceeds 30 microV (threshold +- %d)\r\n', param.threshold));
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'eeglab'
            fileID = fopen(filename, 'a');
            fprintf(fileID, 'data were exported to .set format, following processing steps were performed in EEGLAB\r\n');
            fprintf(fileID, '		- M1 single: dataset with no baseline activity and no zero MEPs --> motor_zero\r\n');
            fprintf(fileID, '		- M1 paired: dataset with no baseline activity --> motor_bl\r\n');
            fprintf(fileID, '		- CTRL: dataset with no post-stimulus motor activity --> motor_zero\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'subset'
            a = find(strcmpi(varargin, 'interpolated'));
            interpolated = varargin{a + 1};
            b = find(strcmpi(varargin, 'discarded'));
            discarded = varargin{b + 1};
            fileID = fopen(filename, 'a');
            fprintf(fileID, '6	data were visually inspected to identify bad channels and epochs\r\n');
            fprintf(fileID, sprintf('       - interpolated channels: %s\r\n', interpolated{1}));
            fprintf(fileID, '       - discarded epochs:\r\n');
            fprintf(fileID, sprintf('               - M1 single:    %d epochs kept\r\n', discarded.epochs_kept(1)));
            fprintf(fileID, sprintf('                               epochs removed: %s\r\n', join(string(discarded.epochs_rejected{1}), ' ')));
            fprintf(fileID, sprintf('               - M1 paired:    %d epochs kept\r\n', discarded.epochs_kept(2)));
            fprintf(fileID, sprintf('                               epochs removed: %s\r\n', join(string(discarded.epochs_rejected{2}), ' ')));
            fprintf(fileID, sprintf('               - CTRL:         %d epochs kept\r\n', discarded.epochs_kept(3)));
            fprintf(fileID, sprintf('                               epochs removed: %s\r\n', join(string(discarded.epochs_rejected{3}), ' ')));            
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'ssp-sir'
            a = find(strcmpi(varargin, 'PC'));
            PC = varargin{a + 1};
            fileID = fopen(filename, 'a');
            fprintf(fileID, '7	muscular artifact was removed using SSP-SIR\r\n');
            fprintf(fileID, '		- first, the filter was calculated and applied to each dataset separately:\r\n');
            fprintf(fileID, sprintf('           - M1 single: %d components removed\r\n', PC(1)));
            fprintf(fileID, sprintf('           - M1 paired: %d components removed\r\n', PC(2)));
            fprintf(fileID, sprintf('           - CTRL: %d components removed\r\n', PC(3)));
            fprintf(fileID, '		- datasets were then additionally filtered by other two filters, to ensure comparability\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
    end
end
function export_EEGLAB(lwdata, filename, ref, subj)
    % dataset
    EEG.setname = filename;
    EEG.filename = [];
    EEG.filepath = [];
    EEG.subject = subj; 
    EEG.session = 1;
    
    % time properties
    EEG.nbchan = lwdata.header.datasize(2);
    EEG.trials = lwdata.header.datasize(1);
    EEG.pnts = lwdata.header.datasize(6);
    EEG.srate = 1/lwdata.header.xstep;
    EEG.times = lwdata.header.xstart + (0:EEG.pnts-1)*lwdata.header.xstep;
    EEG.xmin = EEG.times(1);
    EEG.xmax = EEG.times(end);
    EEG.data = permute(single(lwdata.data),[2,6,1,3,4,5]);
    EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'SEEG_enabled');
    EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'topo_enabled');
    
    % create events with appropriate latencies
    EEG.event = lwdata.header.events;
    if ~isempty(EEG.event)
        [EEG.event.type] = EEG.event.code;
        for e = 1:length(EEG.event)
            EEG.event(e).latency = (e-1)*EEG.pnts + 2001;
        end
        EEG.event = rmfield(EEG.event,'code');
    end
    
    % add specific fields
    EEG.ref = ref;
    EEG.global_tags = lwdata.header.global_tags;
    
    % create required empty fields
    EEG.icawinv = [];
    EEG.icaweights = [];
    EEG.icasphere = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    save([filename,'.set'], 'EEG');
end
function [EEG_out] = SSP_SIR_trials(EEG_in, L, art_topographies, filt_ker, M)
    EEG_out = EEG_in;
    P = eye(size(EEG_in.data,1)) - art_topographies*art_topographies';

    for i = 1:size(EEG_in.data,3)

        data = EEG_in.data(:,:,i);

        %Suppressing the artifacts:
        data_clean = P*data;

        %Performing SIR for the suppressed data:
        PL = P*L;

        if isempty (M)
            M = rank(data_clean) -  1 ;
        end

        tau_proj = PL*PL';
        [U,S,V] = svd(tau_proj);
        S_inv = zeros(size(S));
        S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
        tau_inv = V*S_inv*U';
        suppr_data_SIR = L*(PL)'*tau_inv*data_clean;

        %Performing SIR for the original data:
        tau_proj = L*L';
        [U,S,V] = svd(tau_proj);
        S_inv = zeros(size(S));
        S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
        tau_inv = V*S_inv*U';
        orig_data_SIR = L*(L)'*tau_inv*data;

        if isempty(filt_ker)
            data_correct = suppr_data_SIR;
        else
            filt_ker_B = repmat(filt_ker,[size(suppr_data_SIR,1),1]);
            data_correct = filt_ker_B.*suppr_data_SIR + orig_data_SIR - filt_ker_B.*orig_data_SIR;
            filt_ker_B = [];
        end

        EEG_out.data(:,:,i) = data_correct;

    end
end
