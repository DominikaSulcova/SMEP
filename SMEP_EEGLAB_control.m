%% control script - EEGLAB
% load data, if necessary
addpath(fullfile(folder_toolbox,'eeglab2022.1'));
addpath(fullfile(folder_toolbox,'FastICA_25'));
eeglab
for c = 1:length(condition)
    % load the dataset
    name = sprintf('%s %s S%s visual.set', measure, condition{c}, subj);  
    EEG = pop_loadset('filename', name, 'filepath', folder_data);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    
    % modify dataset info
    EEG.subject = subj; 
    EEG.session = 1;
end
eeglab redraw
%% SSP-SIR
% ACCORDING TO THOMAS:
% determine artifact topograhies from individual datasets
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
    
    % calculate SSP-SIR, use spherical model 
    name = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix); 
    [EEG_temp] = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', time_range, 'PC', []);
    
    % save filter parameters 
    EEG.SSPSIR = EEG_temp.SSPSIR;
    eeglab redraw
    param(c) = EEG.SSPSIR;
    clear EEG_temp
end

% fill in the output structure and save
SMEP.TEP(subject).sspsir_PC = [param(1).PC, param(2).PC, param(3).PC];
SMEP.TEP(subject).sspsir_params = rmfield(param, 'PC');
save([folder_output '\SMEP.mat'], 'SMEP');

% apply ssp-sir using composed filter
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
    
    % visual check - before SSP-SIR
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180]);
    sgtitle(sprintf('S%s - %s: original data', subj, strrep(condition{c}, '_', ' ')))
    figure_name = sprintf('SMEP %s S%s %s no_filter', measure, subj, condition{c}); 
    savefig([folder_output '\Figures\' figure_name '.fig'])
    saveas(gcf, [folder_output '\Figures\' figure_name '.svg'])
            
    % SSP-SIR 
    [EEG] = SSP_SIR_trials(EEG, param(c).L_ave, rmfield(param, {'PC', 'filter', 'L_ave'}), param(c).filter, []);
        
    % baseline correct and save
    EEG = pop_rmbase(EEG, baseline, []);
    name = sprintf('%s %s S%s %s', measure, condition{c}, subj, suffix); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'setname', name, 'overwrite', 'on', 'gui', 'off'); 
    eeglab redraw
    
    % visual check - after SSP-SIR
    figure; 
    pop_timtopo(EEG, [-100  300], [5 15  30  45  60 100 180]);
    sgtitle(sprintf('S%s - %s: filtered data', subj, strrep(condition{c}, '_', ' ')))
    figure_name = sprintf('SMEP %s S%s %s sspsir', measure, subj, condition{c}); 
    savefig([folder_output '\Figures\' figure_name '.fig'])
    saveas(gcf, [folder_output '\Figures\' figure_name '.svg'])
    
    % save dataset
    name = sprintf('%s %s S%s %s.set', measure, condition{c}, subj, suffix); 
    pop_saveset(EEG, 'filename', name, 'filepath', folder_data);
end
%% ICA
% compute ICA matrix on all datasets 
[ALLEEG EEG_all CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',[1:3] ,'study',0); 
EEG_all = eeg_checkset(EEG_all);
EEG = pop_runica(EEG, 'icatype', 'runica', 'concatcond', 'on', 'options', {'extended',1});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% run ICA separately for each dataset
for c = 1:length(condition)
    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
    
    % apply ICA 
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
%% FREQUENCY FILTERS
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