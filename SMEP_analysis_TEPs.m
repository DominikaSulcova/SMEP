% SMEP - TEP preprocessing & analysis
% EEGLAB

%% parameters

%% load the data

fprintf('%s - %s %%rMT:\n', orientation{o}, intensity{i})

% load the dataset
filename = sprintf('AGSICI P1 S%s %s %s.set', subj, orientation{o}, intensity{i});
EEG = pop_loadset('filename', filename, 'filepath', folder_output);
eeglab redraw   

%% baseline correct - subtraction
EEG = pop_rmbase(EEG, baseline,[]);

% visual check
figure; 
pop_timtopo(EEG, [-100  300], [25  45  75  100 180], 'Original data');

%% bad channels
% save the original channels locations 
pop_saveset(pop_select(EEG, 'trial', 1), 'filename', [EEG_all.setname ' all_channels.set'], 'filepath', folder_output);

% visualize the average response to identify possible bad channels
figure; 
pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);

% remove bad channels 
pop_eegplot(pop_select(EEG, 'time', [-0.1 0.3]), 1, 1, 1);
EEG = pop_select(EEG);

%% remove bad trials
%         % verify if the final number of trials is acceptable --> SNR should be less than 10%
%         n_original = 100;
%         n_accepted = 100 - 15;
%         SNR = (sqrt(n_accepted)-sqrt(n_original))/sqrt(n_original);

% remove bad trials
pop_eegplot(EEG, 1, 1, 1);
pop_saveset(EEG, 'filename', [EEG.setname ' trial_select.set'], 'filepath', folder_output);

%% ICA --> remove ocular artifacts and TMS-independent muscular activity
% determine number of components
if EEG.nbchan < 32
    n_comp = ICA_comp - (32 - EEG.nbchan);
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
pop_saveset(EEG, 'filename', EEG.setname, 'filepath', folder_output);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [25  45  75  100 180], 'Data after ICA');

%% SOUND --> suppress extracranial noise 
% - use spherical model
% - regularization level/lambda value can be adjusted to control the
% amount of cleaning - higher value removes more noise, but increases risk of over-correction.
% identify original file with bad channels
chan_file = ['AGSICI P1 S' subj ' ' orientation{o} ' ' intensity{i} ' all_channels.set'];

% run SOUND
EEG = pop_tesa_sound(EEG, 'lambdaValue', lambda, 'iter', 15, 'leadfieldInFile' ,[],  'leadfieldChansFile', [], ...
    'replaceChans', chan_file, 'multipleConds', []);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [25  45  75  100 180], 'Data after SOUND');
pop_saveset(EEG, 'filename', [EEG.setname ' post_SOUND.set'], 'filepath', folder_output);
%         EEG = pop_loadset('filename', ['AGSICI P1 S' subj ' ' orientation{o} ' ' intensity{i} ' pruned with ICA post_ICA.set'],...
%             'filepath', folder_output);

%% SSP-SIR --> leftover muscular artifact
% - use spherical model 
EEG = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', [0,12], 'PC', []);