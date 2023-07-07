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
%       - interpolate the TMS artifact [-0.005 0.002]s
%       - downsample to 2kHz (1/10)
%       - baseine correct by mean subtraction [-250 0.005]
%       - remove extra event categories
% 
% 2) PARSE PER CONDITION
% 
% 3) EXPORT FOR EEGLAB
% 
% 4) LAUNCH EEGLAB
% 
% 5) LOAD DATA TO EEGLAB
% 
% 6) REMOVE BAD CHANNELS & TRIALS
% 
% 7) ICA
%       - ocular artifacts 
%       - TMS-independent muscular activity
%       - electrode noise


%% PARAMETERS
clear all; clc;

% subject
subject = 4;
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% dataset
measure = 'TEP';
block = [1:15];
condition = {'M1_single', 'M1_paired', 'CTRL'}; 

% load output structure 
load([folder_output '\SMEP.mat']);

% create logfile filename, launch the TEP section
filename = sprintf('%s\\Logfiles\\SMEP S%s %s.txt', folder_output, subj, SMEP.info(subject).date); 
logfile_entry('heading', filename);

% choose relevant directories
folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave + eeglab masterfiles
folder_data = uigetdir(pwd, 'Choose the data folder');              % processed data
folder_output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file

% load the finish sound
% load handel
load gong
soundwave = y; clear y Fs

% visualization
fig_counter = 1;

%% 2) PREPROCESSING
% ----- section input -----
suffix = {'reref' 'ep' 'dc' 'interp' 'ds' 'bl'};
epoch = [-1 2];
interp = [-0.005 0.002];
ds_ratio = 10;
baseline = [-0.25 -0.005];
eventcode = 'Stimulation';
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
    [header, data, ~] = RLW_segmentation(header, data, {{eventcode}}, 'x_start', epoch(1), 'x_duration', epoch(2) - epoch(1));
    
    % remove DC and apply linear detrend
    fprintf('removing DC + detrending...')
    [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);

    % interpolate TMS artifact
    fprintf('interpolating TMS artifact...')
    [header, data, ~] = RLW_suppress_artifact_event(header, data,...
        'xstart', interp(1), 'xend',  interp(2), 'event_code', eventcode, 'interp_method', 'pchip');

    % downsample
    fprintf('downsampling...')
    [header, data, ~] = RLW_downsample(header, data, 'x_downsample_ratio', ds_ratio);

    % baseline correct
    fprintf('subtracting baseline...')
    [header, data, ~] = RLW_baseline(header, data, 'operation', 'subtract', 'xstart', baseline(1), 'xend', baseline(2));

    % save for letswave 
    fprintf('done.\n')
    for s = 1:length(suffix)
        header.name = [suffix{s} ' ' header.name];
    end
    CLW_save([], header, data);
end

% create prefix
for s = 1:length(suffix)
    if s == 1
        prefix = suffix{s};
    else
        prefix = [suffix{s} ' ' prefix];
    end
end

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

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

% create prefix
for s = 1:length(suffix)
    if s == 1
        prefix = suffix{s};
    else
        prefix = [suffix{s} ' ' prefix];
    end
end
clear suffix epoch interp ds_ratio baseline eventcode b s data header lwdata option

%% 3) PARSE PER CONDITION
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% load order of stimulation --> stim_order
load([folder_output '\TMS stimulation protocols\SMEP_' subj '_stim_order.mat'])

% label and cluster by condition
fprintf('parsing per condition...\n')
counter = 1;
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

    % merge subjects & save 
    option = struct('type', 'epoch', 'suffix', condition{c}, 'is_save', 1);
    lwdata = FLW_merge.get_lwdata(lwdataset, option); 
end
fprintf('done.\n')
clear b c e option lwdata folder_input stim_order data2merge counter lwdataset 

%% 4) EXPORT FOR EEGLAB
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export in .set format
fprintf('exporting for EEGLAB...\n')
for c = 1:length(condition)
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

%% 6) LOAD DATA TO EEGLAB
% ----- section input -----
this_condition = 'M1_single';
% ------------------------- 
% load the dataset
filename = sprintf('%s S%s %s.set', this_condition, subj, measure);
EEG = pop_loadset('filename', filename, 'filepath', folder_data);
eeglab redraw   

% visual check
figure; 
pop_timtopo(EEG, [-100  300], [15  30  45  60 100 180], 'Original data');

%% 7) REMOVE BAD CHANNELS & TRIALS
% save the original channels locations 
pop_saveset(pop_select(EEG, 'trial', 1), 'filename', [EEG.setname ' all_channels.set'], 'filepath', folder_data);

% visualize the average response to identify possible bad channels
figure; 
pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);

% remove bad channels 
EEG = pop_select(EEG);

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
    end
end


