%% DATA IMPORT
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
% Following script allows to import all three types of data and save it in
% etswave and fieldtrip formates for further preprocessign.

%% parameters
clear all; clc;

% dataset
study = 'SMEP';
subject = 'S04';
block = [1:15];

% choose relevant directories
folder_lw = uigetdir(pwd, 'Choose the letswave folder');        % letswave masterfiles 
folder_output = uigetdir(pwd, 'Choose the output folder');      % processed data
cd(folder_output)

% add letswave 6 to the top of search path
addpath(genpath([folder_lw '\letswave6-master']));

%% import NeurOne data
% define prefix
prefix = {'TEP', 'SMEP'};

% choose folder with raw NeurOne data
folder_input = uigetdir(pwd, 'Coose the input folder');

% specify EEG block subfolders in the imput folder
prompt = {'Admit following EEG blocks:'};
dlgtitle = 'EEG blocks';
dims = [1 35];
definput = {'[1:15]'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
eval(['folders = ' answer ';']); 
clear answer prompt dlgtitle dims definput     

% import the datasets 
for f = 1:length(folders)
    % display current block
    disp([subject  ' - block ' num2str(block(f))])
    
    % import data
    [header, data] = EEG_import_MEGA(folder_input, folders(f));
    
    % fill in letswave history entry
    header.history(1).configuration.parameters.input_folder  = folder_input;
    header.history(1).configuration.parameters.session_number  = folders(f);
    
    % remove extra 'Out' events
    disp('Reducing number of event categories...')
    for i = 1:length(header.events)
        if ~strcmp(header.events(i).code, 'Out')
            index(i) = true;
        else
            index(i) = false;
        end
    end
    header.events = header.events(index);
    fprintf('%d events found in the dataset.', length(header.events))
    clear index
    
    % select TEP data + ECG
    disp('Splitting the datasets & saving')
    data_all = data; header_all = header;
    labels = {header_all.chanlocs([1:30, 50:52]).labels};
    [header, data, ~] = RLW_arrange_channels(header_all, data_all, labels);
    
    % save dataset for letswave
    header.name = [subject ' ' prefix{1} ' b' num2str(block(f))];
    CLW_save([], header, data);
    
    % select SMEP data + ECG, save for letswave
    labels = {header_all.chanlocs([31:52]).labels};
    [header, data, ~] = RLW_arrange_channels(header_all, data_all, labels);
    header.name = [subject ' ' prefix{2} ' b' num2str(block(f))];
    CLW_save([], header, data);
end
clear prefix folder_input folders i f d labels header data header_all data_all

%% import MEP data
% define prefix
prefix = {'MEP'};

% choose the folder with raw data
folder_input = uigetdir(pwd, 'Coose the input folder');

% import the datasets 
for b = block
    % display current block
    disp([subject  ' - block ' num2str(b)])
    
    % import data
    [header, data] = EMG_import_VHDR([folder_input '\' study ' ' subject ' b' num2str(b) '.vhdr']);
    
    % keep only 's1' events
    disp('Reducing number of event categories...')
    for i = 1:length(header.events)
        if strcmp(header.events(i).code, 's1')
            index(i) = true;
        else
            index(i) = false;
        end
    end
    header.events = header.events(index);
    fprintf('%d events found in the dataset.', length(header.events))
    clear index
    
    % discard duplicates if necessary
    if length(header.events) > 80
        fprintf('Removing duplicate triggers...')
        for a = 1:length(header.events)
            latencies(a) = header.events(a).latency;
        end
        [~, index] = unique(latencies, 'stable');
        header.events = header.events(index);
        fprintf('%d unique events found in the dataset.', length(header.events))        
    end
    clear index latencies
    
    % save for letswave
    header.name = [subject ' ' prefix{1} ' b' num2str(b)];
    CLW_save([], header, data);
end
clear prefix folder_input a b i data header


