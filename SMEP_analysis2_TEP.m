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
% ) EXPORT FOR EEGLAB
%       - saves as .set file
% 
% ) LAUNCH EEGLAB
% 
% ) REMOVE MUSCULAR ARTIFACT - SSP-SIR
% 
% ) ICA
%       - removes artifact not timelocked to the stimulus:
%           - ocular artifacts 
%           - TMS-independent muscular activity
%           - electrode noise


%% PARAMETERS
% clear all; clc;

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
% folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave + eeglab masterfiles
% folder_data = uigetdir(pwd, 'Choose the data folder');              % processed data
% folder_output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file

% load output structure 
load([folder_output '\SMEP.mat']);

% create logfile filename, launch the TEP section
filename = sprintf('%s\\Logfiles\\SMEP S%s %s.txt', folder_output, subj, SMEP.info(subject).date); 
logfile_entry('heading', filename);

% load the finish sound
% load handel
load gong
soundwave = y; clear y Fs

% visualization
fig_counter = 1;

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

%% ) EXPORT FOR EEGLAB
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export in .set format
fprintf('exporting for EEGLAB...\n')
for c = 1:length(condition)
    % load the data
    option = struct('filename', sprintf('%s S%s %s', condition{c}, subj, measure));
    lwdata = FLW_load.get_lwdata(option);
    
    % export
    FLW_export_EEGLAB.get_lwdata('filename', sprintf('%s S%s %s', condition{c}, subj, measure), 'pathname', pwd);
end
fprintf('done.\n')
clear c

%% 5) LAUNCH EEGLAB
% add eeglab to the of search path
addpath(fullfile(folder_toolbox,'eeglab2022.1'));

% add fastica to the of search path
addpath(fullfile(folder_toolbox,'FastICA_25'));

% launch eeglab and generate an empty EEGLAB structure
eeglab 

%% 6) REMOVE MUSCULAR ARTIFACT - SSP-SIR
% ----- section input -----
param.time_range = [-0.006, 0.050];
% ------------------------- 
% identify number of opened figures
h = findobj('type','figure'); 
if length(h) > 0
    fig_counter = length(h) + 1;
end
clear h

% cycle through conditions
for c = 1:length(condition)
    % load the dataset
    filename = sprintf('%s S%s %s.set', condition{c}, subj, measure);
    EEG = pop_loadset('filename', filename, 'filepath', folder_data);
    eeglab redraw   
    
    % remove ECG channels 
    EEG = pop_select(EEG, 'channel', {EEG.chanlocs([1:30]).labels});
    eeglab redraw
    
    % visualize the average response to identify possible bad channels
    figure(fig_counter); 
    pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);
    sgtitle(sprintf('S%s - %s: original data', subj, condition{c}))
    fig_counter = fig_counter + 1;
    
    % remove bad channels if necessary
    
    % visual check - before SSP-SIR
    figure(fig_counter); 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Original data');
    sgtitle(sprintf('S%s - %s: original data', subj, condition{c}))
    fig_counter = fig_counter + 1;
    
    % SSP-SIR - spherical model 
    [EEG] = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', param.time_range, 'PC', []);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', sprintf('%s S%s %s SSP', condition{c}, subj, measure),...
        'saveold', sprintf('%s\\%s S%s %s.set', pwd, condition{c}, subj, measure), 'gui', 'off'); 
    eeglab redraw
    
    % visual check - after SSP-SIR
    figure(fig_counter); 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180], 'Filtered data');
    sgtitle(sprintf('S%s - %s: filtered data', subj, condition{c}))
    fig_counter = fig_counter + 1;
end

clear c filename

%% 7) REMOVE BAD CHANNELS & TRIALS
% save the original channels locations 
pop_saveset(pop_select(EEG, 'trial', 1), 'filename', [EEG.setname ' all_channels.set'], 'filepath', folder_data);

% visualize the average response to identify possible bad channels
figure; 
pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);

% remove bad channels 
EEG = pop_select(EEG, 'channel', EEG.chanlocs([1:30]).labels);
eeglab redraw

% remove bad trials
pop_eegplot(EEG, 1, 1, 1);
pop_saveset(EEG, 'filename', [EEG.setname ' trial_select.set'], 'filepath', folder_output);

%% 8) ICA 
% ----- section input -----
ICA_comp = 28;
baseline = [-0.25 -0.005];
% ------------------------- 
% determine number of components
if EEG.nbchan < 30
    n_comp = ICA_comp - (30 - EEG.nbchan);
else
    n_comp = ICA_comp;
end
EEG = pop_tesa_pcacompress(EEG, 'compVal', n_comp, 'plot', 'on');

% run ICA 
EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off');
EEG = pop_tesa_compplot(EEG,'figSize', 'large', 'plotTimeX', [-0.5 0.5], 'plotFreqX', [1 100],...
    'freqScale', 'log', 'saveWeights','off');

% baseline correct & save
EEG = pop_rmbase(EEG, baseline, []);
pop_saveset(EEG, 'filename', EEG.setname, 'filepath', folder_data);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [15  30  45  60 100 180], 'Data after ICA');

%% 9) SOUND --> suppress extracranial noise 
% - use spherical model
% - regularization level/lambda value can be adjusted to control the
% amount of cleaning - higher value removes more noise, but increases risk of over-correction.
% identify original file with bad channels
% ----- section input -----
lambda = 0.025;
% ------------------------- 
chan_file = sprintf('%s S%s %s all_channels.set', this_condition, subj, measure);

% run SOUND
EEG = pop_tesa_sound(EEG, 'lambdaValue', lambda, 'iter', 15);
pop_saveset(EEG, 'filename', [EEG.setname 'post_SOUND.set'], 'filepath', folder_data);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [15  30  45  60 100 180], 'Data after SOUND');

%% 10) SSP-SIR --> leftover muscular artifact
% use spherical model 
EEG = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', [0,12], 'PC', []);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [15  30  45  60 100 180], 'Data after SSP-SIR');

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
    end
end


