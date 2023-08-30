%% SMEP - SCRIPT 3: SMEP preprocessing & analysis
% written by Dominika for the SMEP project (2023)
% 
% The SMEP project explores the recording of spinal motor evoked potentials
% (SMEP) and their use to chanracterize functional changes along the motor
% pathway. To this end, primary motor cortex was stimulated using
% supra-threshold TMS (MagVenture; 120%rMT) and the evoked activity was
% registered at different levels of the cortico-motor pathway:
%   1) TMS-evoked potentials were recorded using 32-channel EEG system
%   (NeurOne, Bittium + EasyCap mounted with TMS compatible multitrodes)
%   2) spinal activity was recorded with a custom 13-channel spinal grid 
%   (multitrodes) + 4 side channels placed on trapezius muscle
%   3) motor-evoked potentials were recorded form the FDI muscle
% 
% Following script provides step-by-step pipeline starting from imported
% raw data (letswave format):
%   1) assign electrode coordinates + detrend
%   ...visual inspection...
%   2) interpolate TMS artifact, interpolate bad electrodes if necessary + remove ECG
%   3) frequency filter, segmentation + linear detrend 
%   4) split into conditions
%   5) visual inspection
%   6) ICA - compute matrix
%   
% 
%% parameters
clear all; clc;

% dataset
measure = 'SMEP';
subject = 'S04';
block = [1:15];
condition = {'M1_single', 'M1_paired', 'CTRL'}; 

% choose relevant directories
folder_git = uigetdir(pwd, 'Choose the Git folder');            % Git repo --> source of saved default files, if needed
folder_toolbox = uigetdir(pwd, 'Choose the letswave folder');   % toolboxes masterfiles
folder_data = uigetdir(pwd, 'Choose the output folder');        % processed data
folder_figures = uigetdir(pwd, 'Choose the figures folder');    % figures

% launch directory and letswave
cd(folder_data)
addpath(genpath([folder_toolbox '\letswave7-master']));
letswave

% visualization
fig_counter = 1;

%% 1) assign electrode coordinates + detrend
% Verify beforehand that you dispose of the file with electrode
% coordinates. If not, disable the first preprocessing step.
% ----- section input -----
suffix = {'dc1'};
coordinates_path = 'F:\letswave7-master\letswave7-master\res\electrodes\spherical_locations\SPINAL_17+1chan.locs';
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% cycle through blocks
for b = block
    % get the dataset
    disp([subject  ' - block ' num2str(b)])
    dataset = [folder_data '\' subject ' ' measure ' b' num2str(b) '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % assign electrode coordinates 
    option = struct('filepath', coordinates_path, 'suffix', '', 'is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    
    % remove DC + linear detrend
    disp('Removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', suffix{1}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
end

% create prefix
prefix = suffix{1};
clear suffix coordinates_path b dataset option lwdata  

%% 2) interpolate TMS artifact, interpolate bad electrodes if necessary + remove ECG
% ----- section input -----
suffix = {'interp' 'ds' 'QRS' 'avg' 'ECG_clean'};
fs = 4000;
duration = 8;
ds_ratio = 5;
interp = [-0.005 0.002];
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% ask for channles to interpolate
chan2interp = inputdlg('Channel to interpolate: ', 'Bad channel', [1 50], {''});
% if length(chan2interp{1}) > 0
%     % identify channels to average
%     keepgoing = 1;
%     chan2avg = cell(0, 0);
%     while keepgoing
%         chan = inputdlg(sprintf('Channel to average #%d: ', length(chan2avg) + 1), 'Bad channel', [1 50], {''});
%         if length(chan{1}) > 0                
%             chan2avg = cat(1, chan2avg, chan);
%         else
%            keepgoing = 0; 
%         end
%     end
% end

% cycle through blocks
for b = block
    disp([subject ' - ' num2str(b)])
    
    % load data and header
    [header, data] = CLW_load([folder_data '\' prefix ' ' subject ' ' measure ' b' num2str(b)]);
    
    % interpolate TMS artifact
    disp('Interpolating the TMS artifact from -0.005 to 0.002 s...')
    [header, data, ~] = RLW_suppress_artifact_event(header, data,...
        'xstart', interp(1), 'xend', interp(2), 'event_code', 'Stimulation', 'interp_method', 'pchip');
    
    % if necessary, interpolate a bad channel
    if length(chan2interp{1}) > 0
        % identify channels to average
        [header, data, ~] = RLW_interpolate_channel(header, data, chan2interp(1), chan2avg);
    end
    
    % downsample
    disp('Downsampling...')
    [header, data, ~] = RLW_downsample(header, data, 'x_downsample_ratio', ds_ratio);
    
    % rereference to right arm
    disp('Re-referencing to right arm...')
    [header, data, ~] = RLW_rereference(header, data, 'apply_list', {'ECGdown'}, 'reference_list', {'ECGright'});
    header_all = header; data_all = data;
    
    % identify QRS in the signal from the rib electrode (II. derivation)
    disp('Identifying QRS complex...')
    [header, data, ~] = RLW_pan_tompkin(header, data, 'channel_label', {'ECGdown'});
    
    % epoch around QRS
    disp('Epoching around QRS...')
    [header, data, ~] = RLW_segmentation(header, data, {'QRS'}, 'x_start', -0.25, 'x_duration', 0.75);
    
    % average & save for letswave
    disp('Saving for letswave...')
    [header, data, ~] = RLW_average_epochs(header, data);
    header.name = [suffix{4} ' ' suffix{3} ' ECG ' prefix ' ' subject ' ' measure ' b' num2str(b)];
    CLW_save([], header, data);

    % rename & clear
    data_ECG = data(1, 1:end - 1, 1, 1, 1, :); 
    header.chanlocs = header.chanlocs(1:end - 1);
    header.datasize = size(data_ECG);
    header_ECG = header;
    clear data header
    
    % compute PCA matrix based on average QRS 
    disp('Computing PCA matrixe...')
    matrix = RLW_PCA_compute(header_ECG, data_ECG);
    
    % visualize components
    data_visual = squeeze(data_all);
    input = matrix.ica_um * data_visual;
    fig = figure(fig_counter);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
    for i = 1:size(data_ECG, 2)
        subplot(4,7,i);
        plot(input(i,(1:(fs * duration))));
        text(3500, 100, num2str(i), 'FontSize', 14, 'FontWeight', 'bold')
        ylim([-150 150])
    end
    fig_counter = fig_counter + 1;

    % decide the number of last components to be removed
    comp2remove = str2num(cell2mat(inputdlg('Number of components to remove: ', 'ECG artifact removal', [1 50], {'[]'})));
    comp2keep = 1:size(data_all, 2);
    comp2keep = comp2keep(~ismember(comp2keep, comp2remove));
    
    % apply PCA matrix, withdraw artifactual components 
    [header_all, data_all] = RLW_ICA_apply_filter(header_all, data_all, matrix.ica_mm, matrix.ica_um, comp2keep);
    
    % save for letswave 
    header_all.name = [suffix{5} ' ' suffix{2} ' ' suffix{1} ' ' header_all.name];
    CLW_save([],header_all, data_all);
end

% update prefix
prefix = [suffix{5} ' ' suffix{2} ' ' suffix{1} ' ' prefix];
clear fs duration ds_ratio interp suffix b i a index apply_list header_all data_all data_ECG header_ECG matrix data_visual input ...
    chan2interp chan2avg keepgoing comp2remove comp2keep fig

%% 3) segmentation   
% ----- section input -----
suffix = {'spinal' 'ep' 'dc2' 'bc1'};
epoch = [-1 1];
baseline = [-0.25 -0.01];
prefix = 'ECG_clean ds interp dc1';
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% process
for b = block
    % get the dataset
    fprintf('%s - block %d: ', subject, b)
    dataset = [folder_data '\' prefix ' ' subject ' ' measure ' b' num2str(b) '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % select spinal electrodes
    fprintf('selecting spinal electrodes ...')
    option = struct('type', 'channel', 'items', {{lwdata.header.chanlocs(2:18).labels}}, 'suffix', suffix{1}, 'is_save', 0);
    lwdata = FLW_selection.get_lwdata(lwdata, option);
    
    % segment
    fprintf('epoching relative to the event ...')
    option = struct('event_labels', {{'Stimulation'}}, 'x_start', epoch(1), 'x_end', epoch(2), 'x_duration', epoch(2) - epoch(1), ...
        'suffix', suffix{2}, 'is_save', 0);
    lwdata = FLW_segmentation.get_lwdata(lwdata, option);
    
    % remove DC + linear detrend
    fprintf('removing DC and applying linear detrend ...')
    option = struct('linear_detrend', 1, 'suffix', suffix{3}, 'is_save', 0);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option); 
    
    % baseline correction
    fprintf('correcting for baseline ...')
    option = struct('operation','substract', 'xstart', baseline(1), 'xend', baseline(2), 'suffix', suffix{4}, 'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata, option);
    fprintf('done.\n')
end

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix epoch baseline b s dataset option lwdata    

%% 4) visual inspection
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));
letswave

% load the output structure
load('SMEP_visual.mat')

% encode manually
visual(str2num(subject(3))).discarded{1} = [];
save('SMEP_visual.mat', 'visual');

%% 5) preliminary ICA
% ----- section input -----
suffix = {'visual' 'prea' 'sp_filter'};
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% load the dataset
for b = block
    dataset{b} = [folder_data '\' suffix{1} ' ' prefix ' ' subject ' ' measure ' b' num2str(b) '.lw6'];
end
option = struct('filename', {dataset});
lwdataset = FLW_load.get_lwdataset(option);

% compute the ICA matrix
fprintf('computing ICA matrix:\n')
fprintf('\n')
option = struct('ICA_mode', 3, 'algorithm', 1, 'percentage_PICA', 100, 'criterion_PICA', 'LAP', 'suffix', suffix{2}, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n')

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix epoch b s dataset option lwdataset 

%% 6) frequency filters 
% ----- section input -----
suffix = {'bandpass' 'notch' 'crop' 'bc2'};
bandpass = [100 1500];
crop = [-0.25 0.75];
baseline = [-0.25 -0.01];
% -------------------------  
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% cycle through blocks
for b = block
    % get the dataset
    fprintf('%s - block %d: ', subject, b)
    dataset = [folder_data '\' prefix ' ' subject ' ' measure ' b' num2str(b) '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % bandpass
    fprintf('applying Butterworth bandpass filter ...')
    option = struct('filter_type', 'bandpass', 'high_cutoff', bandpass(2),'low_cutoff', bandpass(1),...
        'filter_order', 4, 'suffix', suffix{1}, 'is_save', 0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
    
    % notch
    fprintf('applying FFT notch filter ...')
    option = struct('filter_type', 'notch', 'notch_fre', 50, 'notch_width', 2, 'slope_width', 2,...
        'harmonic_num', 5,'suffix', suffix{2},'is_save', 0);
    lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);   
    
    % crop
    fprintf('cropping ...')
    option = struct('xcrop_chk', 1, 'xstart', crop(1), 'xend', crop(2), 'suffix', suffix{3}, 'is_save', 0);
    lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
    
    % baseline correction
    fprintf('correcting for baseline ...')
    option = struct('operation','substract', 'xstart', baseline(1), 'xend', baseline(2), 'suffix', suffix{4}, 'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata, option);
    fprintf('done.\n')
end

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix bandpass crop ds_ratio b s dataset option lwdata 

%% 7) split into conditions
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% load order of stimulation --> stim_order
folder_input = uigetdir(folder_figures, 'Coose the input folder');
load([folder_input '\SMEP_' subject([2:3]) '_stim_order.mat'])

% label and cluster by condition
fprintf('parsing per condition...\n')
counter = 1;
for b = block
    % get the dataset
    dataset = [folder_data '\' prefix ' ' subject ' ' measure ' b' num2str(b) '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % verify the condition
    for c = 1:length(condition)
        if strcmp(condition{c}, stim_order{b})
            % replace event codes
            for e = 1:length(lwdata.header.events)
                lwdata.header.events(e).code = condition{c};
            end
            
            % rename & save into a dataset for merging
            lwdata.header.name = sprintf('%s %s', subject, measure);
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
clear b c e option lwdata folder_input stim_order data2merge counter lwdataset dataset

%% 8) visualize
% set the defaults
channel = 'spinal8';
colours = {[1 0 0] [0.075 0.624 1] [0 0 0]};
int = [-5 2];  
x_lim = [-15, 100];
x_step = 0.5;

% prepare x
x = x_lim(1):x_step:x_lim(2);

% plot 
for c = 1:length(condition)
    % load the data
    load(sprintf('%s S09 SMEP', condition{c}))
end

% launch the figure
fig = figure(figure_counter);
set(gcf, 'units','centimeters','position',[10 10 20 10], 'color', 'w');
hold on
    
% determine axis properties
xl = [x(1) - 25, x(end) + 25];
xlim(xl);
xlabel('time (ms)');  
ylim(y_lim);
ylabel('amplitude (\muV)')

    % loop through channels to plot
    for a = 1:size(data_visual, 1)     
        P(a) = plot(x, data_visual(a, :), 'Color', colour, 'LineWidth', 1);
    end

    % highlight channel
    if ~isempty(channel_n)
        P(end+1) =  plot(x, data_visual(channel_n, :), 'Color', [0 0 0], 'LineWidth', 3);
    end

    % shade interpolated interval 
    rectangle('Position', [int(1), y_lim(1), int(2) - int(1), y_lim(2) - y_lim(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')

    % TMS stimulus
    line([0, 0], y_lim, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % other parameters
    xlim([x(1) - length(x)*(x(2) - x(1))*0.05, x(end) + length(x)*(x(2) - x(1))*0.05])
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
   
%% ICA - compute matrix
% ----- section input -----
suffix = {'ica_PICA'};
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% load the dataset
for c = 1:length(condition)
    dataset{c} = [folder_data '\' prefix2 ' ' condition{c} ' ' subject ' ' measure '.lw6'];
end
option = struct('filename', {dataset});
lwdataset = FLW_load.get_lwdataset(option);

% compute the ICA matrix
option = struct('ICA_mode', 3, 'algorithm', 1, 'percentage_PICA', 100, 'criterion_PICA', 'LAP', 'suffix', suffix, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);

% update prefix 2
prefix2 = [suffix{1} ' ' prefix2];
clear suffix c dataset option lwdataset

%% ICA - extract component features
% ----- section input -----
suffix = {'unmix' 'timecourse' 'frequency'};
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% unmix 
for c = 1:length(condition)
    % get the dataset
    dataset = [folder_data '\' prefix2 ' ' condition{c} ' ' subject ' ' measure '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % unmix the components
    option = struct('suffix', suffix{1}, 'is_save', 0);
    lwdata = FLW_spatial_filter_unmix.get_lwdata(lwdata, option);
    
    % append to the dataset for merging
    lwdataset(c) = lwdata;
end

% merge epochs
disp('Merging epochs...')
option = struct('type', 'epoch', 'suffix', '', 'is_save', 0);
lwdata_merged = FLW_merge.get_lwdata(lwdataset, option);

% average --> timecourse
disp('Averaging...')
option = struct('operation', 'average', 'suffix', suffix{2}, 'is_save', 1);
lwdata = FLW_average_epochs.get_lwdata(lwdata_merged ,option);

% FFT
option = struct('output', 'power', 'half_spectrum', 1, 'suffix', '', 'is_save', 0);
lwdata = FLW_FFT.get_lwdata(lwdata_merged, option);

% average --> frequency content
disp('Averaging...')
option = struct('operation', 'average', 'suffix', suffix{3}, 'is_save', 1);
lwdata = FLW_average_epochs.get_lwdata(lwdata ,option);

% % update prefix
% prefix = ['sp_filter ' prefix];
clear suffix c dataset option lwdataset lwdata lwdata_merged

%% average across blocks
% ----- section input -----
suffix = {'merge_epoch' 'avg_selected'};
% -------------------------  
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% get the datasets
for b = 1:length(block)
    dataset{b} = [folder_data '\' prefix ' ' measure ' ' subject ' ' block{b} '.lw6'];
end
option = struct('filename', {dataset});
lwdataset = FLW_load.get_lwdataset(option);

% merge epochs
disp('Merging epochs...')
option = struct('type', 'epoch', 'suffix', suffix{1}, 'is_save', 0);
lwdata = FLW_merge.get_lwdata(lwdataset, option);

% average epochs
disp('Averaging...')
option = struct('operation', 'average', 'suffix', suffix{2}, 'is_save', 1);
lwdata = FLW_average_epochs.get_lwdata(lwdata ,option);

% update prefix
for s = 1 : length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix b dataset lwdataset lwdata option s

%% calculate correlation
% ----- section input -----
subj = {'pilot_2', 'pilot_4', 'pilot_5'};
time_intervals = {[0.0015, 0.1], [0.0015, 0.015], [0.015, 0.05], [0.05, 0.1]};
time_intervals_labels = {'all', 'early', 'middle', 'late'};
final_dataset = 'avg_selected';
label = 'SELECTED';
% -------------------------  
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% loop through subjects
for s = 1:length(subj)
    % get the single trial data
    data_trials = [];
    for b = 1:length(block)
        dataset = [folder_data '\' prefix ' ' measure ' ' subj{s} ' ' block{b} '.lw6'];
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);
        data_trials = cat(1, data_trials, squeeze(lwdata.data));
    end

    % get the mean data
    dataset = [folder_data '\' final_dataset ' merge_epoch ' prefix ' ' measure ' ' subj{s} ' ' block{1} '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    data_mean = squeeze(lwdata.data);

    % calculate spatial correlation
    for tint = 1:length(time_intervals)
        % select timepoints
        x_start = (time_intervals{tint}(1) - lwdata.header.xstart)/lwdata.header.xstep;
        x_end = (time_intervals{tint}(2) - lwdata.header.xstart)/lwdata.header.xstep;
        x = fix([x_start:x_end]);

        % calculate for each partial average
        for a = 1:500
            if size(data_trials, 1) < a
                % keep the value from previous trial
                corr_spac(s, tint, a) = corr_spac(s, tint, a-1);
            else
                % average trials 
                data_mean_part = squeeze(mean(data_trials(1:a, :, :), 1));

                % calculate spatial correlation at each timepoint
                C = [];
                for t = 1:length(x)
                    u = data_mean(:, x(t));
                    v = data_mean_part(:, x(t));
                    C(t) = (sum(u.*v))/(norm(u)*norm(v));
                end

                % average C across timepoints & append
                corr_spac(s, tint, a) = mean(C);

                % save all timepoints
                if tint == 1
                    corr_spac_all(s, a, :) = C;
                end
            end
        end
    end
end

% define plotting parameters
x = 1:500;
% for s = 1:length(subj)
% 	colours(s, :) = uisetcolor;
% end

% plot spatial correlation by time intervals
for tint = 1:length(time_intervals)
    % launch figure
    fig = figure(fig_counter);
    hold on

    % plot individual data
    for s = 1:length(subj)
        y = squeeze(corr_spac(s, tint, :))';
        p(s) = plot(x, y, 'color', colours(s,:));
        hold on
    end
    
    % plot average across subjects
    corr_spac_mean = squeeze(mean(corr_spac(:, tint, :), 1))';
    p(end + 1) = plot(x, corr_spac_mean, 'color', [0 0 0], 'linewidth', 2.5);
    
    % add cutoff lines
    line([x(1) x(end)], [0.99 0.99], 'linestyle', '--', 'color', [0, 0, 0]) 
    line([x(1) x(end)], [0.95 0.95], 'linestyle', '--', 'color', [0.65, 0.65, 0.65]) 
    
    % calculate cutoff trial numbers and append
    idx_99 = find(corr_spac_mean >= 0.99, 1);
    idx_95 = find(corr_spac_mean >= 0.95, 1);
    
    % mark the cutoff trial numbers
    if ~isempty(idx_99)
        plot(idx_99(1), 0.99, 'marker', 'o', 'markersize', 8, 'markerfacecolor', 'r', 'markeredgecolor', 'none')
        line([idx_99(1) idx_99(1)], [0, 1], 'linestyle', '--', 'color', 'r')
        cutoff(tint, 1) = idx_99(1); 
    end
    if ~isempty(idx_95)
        plot(idx_95(1), 0.95, 'marker', 'o', 'markersize', 8, 'markerfacecolor', 'r', 'markeredgecolor', 'none')
        line([idx_95(1) idx_95(1)], [0, 1], 'linestyle', '--', 'color', 'r')
        cutoff(tint, 1) = idx_95(1); 
    end
    
    % adjust plot visuals
    set(gca, 'fontsize', 12)
    xlabel('averaged trials')
    ylabel('spatial correlation')
    ylim([0 1])
    legend(p, [subj {'mean'}], 'location', 'southwest', 'fontsize', 12, 'edgecolor', 'none', 'numcolumns', 2)
    title(['spatial correlation: ' time_intervals_labels{tint}])
    
    % save figure
    fig_name = ['spat_corr_' label '_' time_intervals_labels{tint}];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'])
    fig_counter = fig_counter + 1;
    
    % clear workspace
    clear p
end

clear subj time_intervals time_intervals_labels b dataset option data_trials s tint x_start x_end x y ...
    a u v t C corr_spac_mean 9 idx_99 idx_95 fig_name p data_visual fig label

%%
% plot image with scaled colors
data_visual = squeeze(mean(corr_spac_all, 1))';
x = 1.5 : lwdata.header.xstep*1000 : 100;
y = 1:size(data_visual, 2);

% launch figure
fig = figure(fig_counter);
hold on
imagesc(x, y, data_visual)
colorbar

% clear data_visual x y fig


%% plot scMEP timecourse and peak topography
% ----- section input -----
time_window = [0 0.025];
shading = [0.002 0.012];
% -------------------------
% load the data & header
load([folder_data '\' prefix ' ' measure ' ' subject ' ' block{1} '.mat']);
load([folder_data '\' prefix ' ' measure ' ' subject ' ' block{1} '.lw6'], '-mat');

% plot timecourse
t = header.xstart + [header.xstep*(1 : header.datasize(6))];
fig = figure(fig_counter);
hold on
plot_timecourse(data, header, time_window, 'shading', shading)
hold off

% save figure
fig_name = [measure ' ' subject ' timecourse'];
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'])
fig_counter = fig_counter + 1;

% identify peak features for topography plots
finish = 0;
counter = 1;
while finish == 0
    % wait until the mouse is clicked
    w = waitforbuttonpress;

    % get the position of the mouse
    CP = get(gca, 'CurrentPoint');
    
    % append it to the peak vector
    peak_latency(counter) = CP(1,1);
    peak_amplitude(counter) = CP(1,2);
    counter = counter + 1;
    
    % ask for approval
    answer = questdlg('Do you want to select another peak?', 'Peak latency selection',...
        'Yes.', 'No, proceed to plotting.', 'No, proceed to plotting.');    
    switch answer
        case 'Yes.'
        case 'No, proceed to plotting.'
            % exit the while loop
            finish = 1;
    end
    
end
clear w CP counter answer

% plot topography
for a = 1:length(peak_latency)
    % identify data
    [~, t_idx] = min(abs(t - peak_latency(a)));
    v = double(squeeze(data(1,1:17,1,1,1,t_idx)));

    % spinal topoplot coordinates
    x =([-1 0 1 -1 0 1 -1 0 1 -1 0 1 0 -2.5 2.5 -2.5 2.5]);
    y =([3 4 3 1 2 1 -1 0 -1 -3 -2 -3 -4 2 2 -2 -2]);
    
    % determine map limits
    map_lim = [-1 * abs(peak_amplitude(a)), abs(peak_amplitude(a))];

    % plot topography
    fig = figure(fig_counter);
    hold on
    plot_topography(x, y, v, 'maplimits', map_lim)
    hold off

    % save figure
    fig_name = [measure ' ' subject ' ' num2str(a) ' topography'];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'])
    fig_counter = fig_counter + 1;
end

clear time_window shading map_lim b data header a t t_idx v x y fig fig_name 

%% muscular artifact
% ----- section input -----
suffix = {'muscle' 'bandpass' 'notch' 'ep' 'dc2' 'avg'};
bandpass = [2 1800];
epoch = [-0.05 0.1];
prefix_old = 'ECG_clean ds interp dc1';
coordinates_path = 'F:\letswave7-master\letswave7-master\res\electrodes\spherical_locations\SPINAL_17chan.locs';
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% prepare the dataset
for b = 1:length(block)
    % get the data
    disp([subject ' - ' block{b}])
    dataset = [folder_data '\' prefix_old ' ' measure ' ' subject ' ' block{b} '.lw6'];
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
    
    % assign electrode coordinates 
    option = struct('filepath', coordinates_path, 'suffix', '','is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    
    % select side electrodes
    disp('Selecting side electrodes...')
    option = struct('type', 'channel', 'items',{electrodes([14:18, 22, 23])}, 'suffix', suffix{1}, 'is_save', 0);
    lwdata = FLW_selection.get_lwdata(lwdata, option);
    
    % bandpass
    disp('Applying Butterworth bandpass filter...')
    option = struct('filter_type', 'bandpass', 'high_cutoff', bandpass(2),'low_cutoff', bandpass(1),...
        'filter_order', 4, 'suffix', suffix{2}, 'is_save', 0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
    
    % notch
    disp('Applying FFT notch filter...')
    option = struct('filter_type', 'notch', 'notch_fre', 50, 'notch_width', 2, 'slope_width', 2,...
        'harmonic_num', 2,'suffix', suffix{3},'is_save', 0);
    lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);
    
    % segment
    disp('Epoching relative to the event...')
    option = struct('event_labels', {{'Stimulation'}}, 'x_start', epoch(1), 'x_end', epoch(2), 'x_duration', epoch(2) - epoch(1), ...
        'suffix', suffix{4}, 'is_save', 0);
    lwdata = FLW_segmentation.get_lwdata(lwdata, option);
    
    % remove DC + linear detrend
    disp('Removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', suffix{5}, 'is_save', 0);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option); 
    
    % append to the dataset for merging
    lwdataset(b) = lwdata;
end

% merge epochs
disp('Merging epochs...')
option = struct('type', 'epoch', 'suffix', '', 'is_save', 0);
lwdata = FLW_merge.get_lwdata(lwdataset, option);

% average epochs
disp('Averaging...')
option = struct('operation', 'average', suffix{2}, 'avg', 'is_save', 1);
lwdata = FLW_average_epochs.get_lwdata(lwdata ,option);

clear suffix prefix_old coordinates_path bandpass epoch b dataset option lwdata lwdataset 

%% functions
function plot_timecourse(data, header, time_window, varargin)   
    % identify time axis and voltage 
    t_visual = time_window(1) : header.xstep : time_window(2);
    t_start = (time_window(1) - header.xstart)/header.xstep;
    t_end = (time_window(2) - header.xstart)/header.xstep;
    v_visual = double(squeeze(data(1,8,1,1,1,t_start : t_end)));
    
    % check for y limits
    if ~isempty(varargin)
        a = find(strcmpi(varargin, 'y_lim'));
        if ~isempty(a)
            y_lim = varargin{a + 1};
        else
            plot(t_visual, v_visual)
            y_lim = get(gca, 'ylim');
            cla
        end
    end    
    
    % check for shading limits
    if ~isempty(varargin)
        b = find(strcmpi(varargin, 'shading'));
        if ~isempty(b)
            shading_lim = varargin{b + 1};
        else
            shading_lim = [];
        end
    end  
    
    % shade MEP latency
    if ~isempty(shading_lim)
        rectangle('position', [shading_lim(1), y_lim(1), shading_lim(2)-shading_lim(1), y_lim(2)-y_lim(1)], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')
    end
    
    % plot the data
    plot(t_visual, v_visual, 'Color', [0 0 0], 'LineWidth', 2.5)
    
    % figure parameters
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    xlim([time_window(1) - 0.001, time_window(2) + 0.001])
    ylim([y_lim])
    set(gca, 'FontSize', 14) 
    set(gca, 'layer', 'top');
    
    % identify the channel
    text(0.015, y_lim(2) - (y_lim(2) - y_lim(1))/10, 'electrode at C6', 'FontSize', 14);
end
function plot_topography(x, y, v, varargin) 
    % determîne map limits
    if ~isempty(varargin)
        a = find(strcmpi(varargin, 'maplimits'));
        if ~isempty(a)
            maplimits = varargin{a + 1};
        else
            maplimits(1) = min(v);
            maplimits(2) = max(v);
        end
    end

    % plot the data
    xi = linspace(-2.5, 2.5, 1000);
    yi = linspace(-4, 4 , 1000);
    [Xi,Yi,Zi] = griddata(x, y, v, xi', yi,'cubic');
    surface(xi, yi, zeros(size(Zi)), Zi, 'EdgeColor', 'none','FaceColor', 'flat');
    plot(x, y, 'k.', 'MarkerSize', 25);
    contour(Xi, Yi, Zi, 5, 'k');
    
    % figure parameters
    axis off;
    colorbar('southoutside', 'box', 'off');
    set(gcf, 'position', [100,100,300,450]);
    set(gca,'clim',maplimits);

    % highlight C6
    plot(x(8), y(8), 'r.', 'MarkerSize', 25)
    text(x(8) + 0.15, y(8), 'C6');
end