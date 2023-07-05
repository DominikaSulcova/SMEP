%% SMEP - SCRIPT 1: MEP preprocessing & analysis
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
% Following script provides the code for complete processing of MEPs
% - each section of the script performs one processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - output variables are gathered in an output structure 'SMEP.mat'
% 
% 1) SUBJECT & SESSION INFORMATION 
%       - launch info structure 'SMEP.mat' and creates a field 'info'
%       - encodes subject characteristics: age, sex, handedness, rMT,
%       start/end of the experimental session, stimulation sites (Visor2)
%       - initiates the logfile, fills in the info
%
% 2) PREPROCESSING
%       - segmentation relative to events --> epochs [-0.3 0.3]s, epoch size 614 bins
%       - DC removal + linear detrend
% 
% 3) REMOVE MISSED OR FAULTY TRIALS
%       - requires manual input of missed trials: a cell 'missed' 
%           - first column indicates block number
%           - second column contains a vector of indexes of missed trials
% 
% 4) PARSE PER CONDITION
%       - loads the stimulation order (needs to be in the data folder!)
%       - labels the epochs with appropriate eventcodes
%       - merges epochs into three condition datasets 
% 
% 5) DISCARD TRIALS WITH EXCESSIVE BASELINE ACTIVITY

%% PARAMETERS
clear all; clc;

% dataset
measure = 'MEP';
subject = 4;
block = [1:15];
condition = {'M1_single', 'M1_paired', 'CTRL'}; 

% identify subject
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% choose relevant directories
folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
folder_data = uigetdir(pwd, 'Choose the data folder');              % processed data
folder_output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file

% load the finish sound
% load handel
load gong
soundwave = y; clear y Fs

% visualization
fig_counter = 1;

%% 1) SUBJECT & SESSION INFORMATION
% load if available
if exist([folder_output '\SMEP.mat']) > 0 
    load([folder_output '\SMEP.mat']);
end

% ask for subject information
prompt = {'Age:', 'Sex:', 'Handedness:'};
dlgtitle = sprintf('Subject %d: personal information', subject);
dims = [1 25];
definput = {'25', 'female', 'right'};
info = inputdlg(prompt, dlgtitle, dims, definput);
clear prompt dlgtitle dims definput

% encode subject information to the output structure
SMEP.info(subject).ID = subject;
SMEP.info(subject).age = str2num(info{1});
SMEP.info(subject).sex = info{2};
SMEP.info(subject).handedness = info{3};

% ask for session information
prompt = {'Date:', 'Start:', 'End:', 'rMT:', 'Motor hotspot:', 'Control site:'};
dlgtitle = sprintf('Subject %d: session information', subject);
dims = [1 25];
definput = {'00002023', '9:00', '13:00', '50', '00', '00'};
info = inputdlg(prompt, dlgtitle, dims, definput);
clear prompt dlgtitle dims definput

% encode session information to the output structure
SMEP.info(subject).date = info{1};
SMEP.info(subject).start = info{2};
SMEP.info(subject).end = info{3};
SMEP.info(subject).rMT = str2num(info{4});
SMEP.info(subject).M1 = str2num(info{5});
SMEP.info(subject).CTRL = str2num(info{6});

% save the output structure
save([folder_output '\SMEP.mat'], 'SMEP');

% create logfile
filename = sprintf('%s\\Logfiles\\SMEP S%s %s.txt', folder_output, subj, SMEP.info(subject).date); 
initialize_logfile(SMEP, subject, filename);
logfile_entry('heading', filename);
clear info 

%% 2) PREPROCESSING
% ----- section input -----
suffix = {'ep' 'dc'};
epoch = [-0.3 0.3];
eventcode = 's1';
% -------------------------
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% preprocess in letswave
for b = block
    fprintf('block %d:\n', b)

    % load data and header
    [header, data] = CLW_load(sprintf('S%s %s b%d', subj, measure, b));   
    
    % segment 
    fprintf('epoching...')
    [header, data, ~] = RLW_segmentation(header, data, {{eventcode}}, 'x_start', epoch(1), 'x_duration', epoch(2) - epoch(1));
    
    % remove DC and apply linear detrend
    fprintf('removing DC + detrending...')
    [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);

    % save for letswave 
    fprintf('done.\n')
    for s = 1:length(suffix)
        header.name = [suffix{s} ' ' header.name];
    end
    CLW_save([], header, data);
end

create prefix
for s = 1:length(suffix)
    if s == 1
        prefix = suffix{s};
    else
        prefix = [suffix{s} ' ' prefix];
    end
end

% update logfile
logfile_entry('preprocessing', filename);
clear suffix epoch eventcode b s data header

%% 3) REMOVE MISSED OR FAULTY TRIALS
% ----- section input -----
missed = {1, [1];...
        9, [7]; ...
        11, [1:80]; ... 
        12, [1:29]}; 
% -------------------------
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% append the missed epochs to the output file
SMEP.MEP(subject).ID = subject;
SMEP.MEP(subject).missed_epochs = missed;

% discard missing/faulty epochs if necessary
if ~isempty(missed)
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

        % discard indicated epochs 
        option = struct('type', 'epoch', 'items', {epochs2keep}, 'suffix', '', 'is_save', 1);
        lwdata = FLW_selection.get_lwdata(lwdata, option);    
    end
    
    % encode to the logfile 
    logfile_entry('missed', filename, 'missed_epochs', missed);
else    
    % encode to the logfile 
    logfile_entry('missednot', filename);
end

% save the output structure
save([folder_output '\SMEP.mat'], 'SMEP');
clear missed b option lwdata epochs epochs2keep e

%% 4) PARSE PER CONDITION
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
    
%     % in case a whole block needs to be removed
%     if c == 1
%        lwdataset = lwdataset([1:3, 5]);
%     end

    % merge subjects & save 
    option = struct('type', 'epoch', 'suffix', condition{c}, 'is_save', 1);
    lwdata = FLW_merge.get_lwdata(lwdataset, option); 
end
fprintf('done.\n')
clear b c e option lwdata folder_input stim_order data2merge counter lwdataset 

%% 5) DISCARD TRIALS WITH EXCESSIVE BASELINE ACTIVITY
% ----- section input -----
prefix = 'visual';
baseline = [-0.2 0];
threshold = 15;
allow_sd = 3;
% -------------------------
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% loop through conditions
for c = 1:length(condition)                              
    fprintf('condition: %s\n', condition{c})
    
    % load data and header
    [header, data] = CLW_load(sprintf('%s S%s %s', condition{c}, subj, measure));   

    % prepare variables 
    x = baseline(1) : header.xstep : 0; 
    data_visual = squeeze(data(:, 1, 1, 1, 1, ceil((baseline(1) - header.xstart) / header.xstep) : floor((baseline(2) - header.xstart) / header.xstep)));
    data_visual_orig = data_visual;
    discarded = [];

    % ------- run automatic RMS + SD removal ------- 
    go = 1; cycle = 1; 
    while go 
        discarded_pos = [];

        % calculate average RMS values
        avg_rms = mean(rms(data_visual'));
        avg_sd = std(rms(data_visual'));
        cutoff = avg_rms + allow_sd * avg_sd;

        % plot original dataset
        fig = figure(fig_counter); 
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        subplot(3, 2, [1 2])
        plot(x, data_visual, 'color', [0.45, 0.45, 0.45], 'linewidth', 1.5)                
        sgtitle(sprintf('Subject %d: CYCLE %d', subject, cycle), 'fontsize', 18)
        title(sprintf('Original dataset - average RMS %d', avg_rms))
        set(gca,'Fontsize',16); ylabel('amplitude (\muV)');
        hold on

        % loop through trials
        for e = 1:length(header.events)
            % check if the rms of current event fits into the limis 
                current_rms = rms(data_visual(e, :)');
                if current_rms > cutoff
                    discarded_pos = [discarded_pos e];
                    discarded = [discarded header.events(e).epoch];
                    ditch = true;
                else
                    ditch = false;
                end                        

            % plot the event to the appropriate axes                
            if ditch   
                subplot(3, 2, [5 6])
                plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
                hold on
            else
                subplot(3, 2, [3 4])
                plot(x, data_visual(e, :), 'Color', [0, 0, 0], 'LineWidth', 1.5)
                hold on
            end
        end

        % add parameters to axes
        subplot(3, 2, [3 4])
        title(['Kept epochs: ' num2str(length(header.events) - length(discarded_pos))])
        set(gca,'Fontsize',16); ylabel('amplitude (\muV)');
        hold on

        subplot(3, 2, [5 6])
        title(['Discarded epochs: ' num2str(length(discarded_pos))])
        set(gca,'Fontsize',16); xlabel('time (s)'); ylabel('amplitude (\muV)');
        hold on

        % cycle to 0 discarded epochs 
        if ~isempty(discarded_pos)
            % remove indicated epochs from data, update header
            data(discarded_pos, :, :, :, :, :) = [];
            header.datasize(1) = header.datasize(1) - length(discarded_pos);
            header.events(discarded_pos)= [];

            % remove indicated epochs from visual dataset for
            % future filtration cycles
            data_visual(discarded_pos , :) = [];

            % continue the cycle
            pause(1); clf;
            cycle = cycle + 1;  
        else
            % close current figure
            pause(1); close(fig);   

            % plot original dataset
            fig = figure(fig_counter); 
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            subplot(3, 2, 3)
            plot(x, data_visual_orig, 'color', [0, 0, 0], 'linewidth', 1.5)                       
            title(['Original dataset: ' num2str(size(data_visual_orig, 1)) ' epochs'])
            set(gca,'Fontsize',14); ylabel('amplitude (\muV)');
            hold on  

            % plot filtered dataset
            subplot(3, 2, 5)
            plot(x, data_visual, 'Color', [0, 0, 0], 'linewidth', 1.5)
            xlim = get(gca,'xlim');                        
            title(['RMS + SD: ' num2str(length(header.events)) ' epochs kept, ' num2str(cycle - 1) ' cycles performed'])
            set(gca,'Fontsize',14); ylabel('amplitude (\muV)'); xlabel('time (s)');
            hold on

            % add threshold 
            subplot(3, 2, 5)
            l(1) = line([-0.2, 0], [threshold, threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
            l(2) = line([-0.2, 0], [-threshold, -threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
            text(xlim(1) + 0.005 , - threshold + 4 ,['threshold = ' num2str(threshold) ' \muV'], 'Fontsize', 14, 'color', [0.99, 0.3, 0.2])
            hold on

            % plot discarded epochs - RMS + SD
            if ~isempty(discarded)
                subplot(3, 2, 4)
                plot(x, data_visual_orig(discarded, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5) 
                title(['RMS + SD: ' num2str(length(discarded)) ' epochs discarded'])
                set(gca,'Fontsize',14)
                hold on
            else
                subplot(3, 2, 4)
                title('No epoch discarded'); set(gca,'Fontsize',14)
                hold on
            end                     

            % exit the while loop
            go = 0;
        end
    end

    % ------- remove epochs that depass the threshold -------
    % loop through left trials
    for e = 1:length(header.events)
        % check if tthe maximum value across baseline datapoints
        % fits under the threshold
            current_max = max(abs(data_visual(e, :)));
            if current_max > threshold
                discarded_pos = [discarded_pos e];
                discarded = [discarded header.events(e).epoch];
                ditch = true;
            else
                ditch = false;
            end                        

        % plot the event to the appropriate axes                
        if ~ditch   
            subplot(3, 2, [1 2])
            plot(x, data_visual(e, :), 'Color', [0, 0.45, 0.74], 'LineWidth', 1.5)
            hold on
        else
            subplot(3, 2, 6)
            plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
            hold on
        end
    end

    % remove indicated epochs from visual dataset
    data_visual(discarded_pos , :) = [];

    % add parameters to axes
    final_rms = mean(rms(data_visual'));
    subplot(3, 2, [1 2])
    set(gca,'Fontsize',14); ylabel('amplitude (\muV)');  
    title(sprintf('Subject %d - FINAL DATASET: %d epochs kept, %d discarded - final average RMS %d', ...
        subject, length(header.events) - length(discarded_pos), length(discarded), final_rms), 'Fontsize', 18)            
    hold on

    subplot(3, 2, 6)
    title(['THRESHOLD : ' num2str(length(discarded_pos)) '  epochs discarded'])
    set(gca,'Fontsize',14); xlabel('time (s)');
    hold on

    pause(3)

    % ------- save outcome variables -------          
    % modify and save header
    header.datasize(1) = header.datasize(1) - length(discarded_pos);
    header.events(discarded_pos) = [];
    header.events = header.events';
    for h = 1:header.datasize(1)
        header.events(h).epoch = h;
    end
    header.name = [prefix ' ' header.name];
    save([header.name '.lw6'], 'header')         

    % modify and save data
    data(discarded_pos, :, :, :, :, :) = [];
    save([header.name '.mat'], 'data')

    % fill in the outcome table and save
    SMEP.MEP(subject).baseline_discarded{c} = sort(discarded);
    SMEP.MEP(subject).baseline_cycles(c) = cycle - 1;
    SMEP.MEP(subject).baseline_threshold = threshold;

    % save the output structure
    save([folder_output '\SMEP.mat'], 'SMEP');

    % save and close the figure
    figure_name = sprintf('SMEP %s baseline S%s %s', measure, subj, condition{c}); 
    savefig([folder_output '\Figures\' figure_name '.fig'])
    saveas(fig, [folder_output '\Figures\' figure_name '.svg'])
    close(fig)

    % update the counter
    fig_counter = fig_counter + 1;
end

%% update the logfile
logfile_entry('baseline', filename)
clear baseline threshold allow_sd c header data data_visual x 

%% FUNCTIONS
function  initialize_logfile(SMEP, subject, filename)
    % identify subject
    if subject < 10
       subj = ['0' num2str(subject)];
    else
       subj = num2str(subject); 
    end

    % fill in the file
    fileID = fopen(filename, 'w');
    fprintf(fileID, '******************************************************************************************************\r\n');
    fprintf(fileID, sprintf('study: SMEP\r\n')); 
    fprintf(fileID, sprintf('subject: S%s\r\n', subj));
    fprintf(fileID, sprintf('         %s\r\n', SMEP.info(subject).sex));
    fprintf(fileID, sprintf('         %d y.o.\r\n', SMEP.info(subject).age));
    fprintf(fileID, sprintf('         %shanded\r\n', SMEP.info(subject).handedness));
    fprintf(fileID, sprintf('date: %s/%s/%s\r\n', SMEP.info(subject).date(1:2), SMEP.info(subject).date(3:4), SMEP.info(subject).date(5:8)));
    fprintf(fileID, '******************************************************************************************************\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, 'DATA ACQUISITION\r\n');
    fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '- TMS: MagVenture MagPro X100 + Visor2 neuronavigation based on default brain model\r\n');
    fprintf(fileID, '     - biphasic sin pulse, PA-AP current direction\r\n');
    fprintf(fileID, '     - 3 tested conditions:    1) single pulse over left M1 - 120 %%rMT\r\n')
    fprintf(fileID, '                               2) SICI over left M1 - 80 + 120 %%rMT, 2.5ms ISI\r\n')
    fprintf(fileID, '                               3) single pulse over control site - 120 %%rMT\r\n')
    fprintf(fileID, '- EEG + SEG: Bittium system, TESLA amplifiers, NeurOne recording software\r\n');
    fprintf(fileID, '     - 20kHz sampling rate, 3500Hz LP device filter\r\n');
    fprintf(fileID, '     - referenced to linked earlobes, ground on the right cheek\r\n');
    fprintf(fileID, '     - EEG --> amplifier 1, 30 electrodes (10-20), M1/M2 not recorded\r\n');
    fprintf(fileID, '     - SEG --> amplifier 2, 19 electrodes: 13 spinal, 4 side, epiglottis, thoracic\r\n');
    fprintf(fileID, '     - ECG --> amplifier 2, 3 electrodes: left arm, left ribs, right arm\r\n');
    fprintf(fileID, '- recorded in 15 blocks:\r\n');
    fprintf(fileID, '     - each block consisted of 80 stimuli of the same condition with jittered ITI 3 - 5s\r\n');
    fprintf(fileID, '     - the sequence of conditions was pseudorandom --> split in shuffled triplets\r\n');
    fprintf(fileID, '     - each triplet of blocks was delivered by different experimentator - Dominika or Alessandra\r\n');
    fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
    fprintf(fileID, '\r\n');
    fprintf(fileID, sprintf('session start: %s\r\n', SMEP.info(subject).start));
    fprintf(fileID, sprintf('session end: %s\r\n', SMEP.info(subject).end));
    fprintf(fileID, '\r\n');
    fprintf(fileID, sprintf('rMT: %d %%rMT\r\n', SMEP.info(subject).rMT));
    fprintf(fileID, '\r\n');
    fprintf(fileID, sprintf('motor hotspot: stim %d\r\n', SMEP.info(subject).M1));
    fprintf(fileID, sprintf('control site: stim %d\r\n', SMEP.info(subject).CTRL));
    fprintf(fileID, '\r\n');
    fprintf(fileID, 'notes:\r\n');
    fprintf(fileID, '\r\n');
    fclose(fileID);
end
function logfile_entry(entry, filename, varargin)
    switch entry
        case 'heading'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, 'MOTOR-EVOKED POTENTIALS\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, '1	dataset imported using matlab script ''SMEP_import.m''\r\n');
            fprintf(fileID, '		- VHDR data loaded from Visor2 output\r\n');
            fprintf(fileID, '		- only ''s1'' category kept\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);

        case 'preprocessing'
            fileID = fopen(filename, 'a');
            fprintf(fileID, 'following processing steps were performed using matlab script ''SMEP_analysis1_MEP.m''\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '2	preprocessing\r\n');
            fprintf(fileID, '		- segmentation relative to events --> epochs [-0.3 0.3]s, epoch size 614 bins\r\n');
            fprintf(fileID, '		- DC removal + linear detrend\r\n');
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

        case 'ICA 1'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'components'));
                if ~isempty(a)
                    components = varargin{a + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '4	preliminary ICA was computed to get rid of the major decay artifact\r\n');
            fprintf(fileID, '		- squared matrix --> 32 ICs\r\n');
            fprintf(fileID, '		- removed component(s): IC %s\r\n', num2str(components));
            fprintf(fileID, '	--> file prefix: sp_filter prea\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 2'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '5	frequency filters applied\r\n');
            fprintf(fileID, '		- FFT notch filter --> 50 Hz, width 2 + slope 2\r\n');
            fprintf(fileID, '		- Butterworth bandpass filter --> [0.1 80]Hz, 4th order\r\n');
            fprintf(fileID, '	--> file prefix: bandpass notch\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '6 	signal cropped to get rid of edge artifacts\r\n');
            fprintf(fileID, '		- [-0.75 0.75]s\r\n');
            fprintf(fileID, '	--> file prefix: crop\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'visual'
            if ~isempty(varargin)
                e = find(strcmpi(varargin, 'epochs'));
                if ~isempty(e)
                    epochs = varargin{e + 1};
                end
                
                n = find(strcmpi(varargin, 'epochs_n'));
                if ~isempty(n)
                    epochs_n = varargin{n + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '7	noisy epochs were removed based on visual inspection\r\n');
            fprintf(fileID, '		- spTMS TS: %s --> %d epochs retained\r\n', num2str(epochs(1, epochs(1, :) > 0)), epochs_n(1));
            fprintf(fileID, '		- spTMS CS1: %s --> %d epochs retained\r\n', num2str(epochs(2, epochs(2, :) > 0)), epochs_n(2));
            fprintf(fileID, '		- spTMS CS2: %s --> %d epochs retained\r\n', num2str(epochs(3, epochs(3, :) > 0)), epochs_n(3));
            fprintf(fileID, '		- spTMS CS3: %s --> %d epochs retained\r\n', num2str(epochs(4, epochs(4, :) > 0)), epochs_n(4));
            fprintf(fileID, '		- ppTMS CS1: %s --> %d epochs retained\r\n', num2str(epochs(5, epochs(5, :) > 0)), epochs_n(5));
            fprintf(fileID, '		- ppTMS CS2: %s --> %d epochs retained\r\n', num2str(epochs(6, epochs(6, :) > 0)), epochs_n(6));
            fprintf(fileID, '		- ppTMS CS3: %s --> %d epochs retained\r\n', num2str(epochs(7, epochs(7, :) > 0)), epochs_n(7));
            fprintf(fileID, '	--> file prefix: visual\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 3'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '9  baseline corrected - subtracted mean of [-0.2 -0.005]s\r\n');
            fprintf(fileID, '   --> file prefix: bl_correct\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '10 signal averaged across epochs\r\n');
            fprintf(fileID, '   --> file prefix: avg\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
    end
end