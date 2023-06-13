%% Stimulation protocol generator for MagPro X100 
% Writtent by Domi for the SMEP project (2023)
% 
% 1)    two stimulation protocols with jittered inter-stimulus interval 
%       (3 - 5s ISI)are generated for each subject:
%           - 80 single-pulse TMS stimuli, supra-threshold intensity - 120%rMT
%           - 80 paired-pulse TMS stimuli - 2.5ms ISI, 80 + 120%rMT
% 
% 2)    a sequence of 15 stimulation blocks was generated to randomize and
%       balance the 3 conditions within each subject and across subjects:
%           - single-pulse over the M1 (M1_single)
%           - paired-pulse over the M1 (M1_paired)
%           - single-pulse over the control site (CTRL)

%% parameters
study = 'SMEP';
subject = 'S02';

% stimulation parameters
amplitude_TS = 120;                 % testing stimulus intensity --> %rMT
amplitude_CS = 80;                  % conditioning stimulus intensity
repetitions = 80;                   % number of trials in a block
waiting_time = [3000 6000];         % interval disponible for ITI
ISI = 2500;                         % inter stimulus interval for SICI

% sequence parameters
subjects = [2:21];
blocks = 15;
conditions = {'M1_single', 'M1_paired', 'CTRL'};

%% write the protocols
% single-pulse, supra-threshold
wait_sequence = (randi([waiting_time(1) waiting_time(2)], repetitions, 1)); 
filename = strcat(subject, '_single.CG3');
initializeMagProFilebi(filename, repetitions); 
for i= 1:repetitions
    if i==1
        delay = 10000;
    else
        delay = wait_sequence(i);
    end
    writeMagPro_singlepulsebi(filename, i, delay, amplitude_TS);
end  

% paired-pulse
wait_sequence = (randi([waiting_time(1) waiting_time(2)],repetitions, 1)); 
filename = strcat(subject, '_paired.CG3');
initializeMagProFilebi(filename,repetitions); 
for i= 1:repetitions
    if i==1
        delay = 10000;
    else
        delay = wait_sequence(i);
    end
    MagPro_paired_bi_R(filename, i, delay, amplitude_CS, amplitude_TS, ISI);
end
clear wait_sequence filename i delay 

%% generate protocol sequence 
% identify all permutations of the conditions
perm_cond = perms(conditions);

% loop through subjects
stim_order_all = {};
for s = subjects
    % add 0 if s < 10
    if s < 10
        subj_name = ['0', num2str(s)];
    else
        subj_name = num2str(s);
    end
    
    % initialize the output matrix 
    stim_order = cell(blocks, 1);
    sequence = randperm(size(perm_cond, 1));
    for b = 1:blocks/length(conditions)
        stim_order((b - 1)*length(conditions) + 1 : (b - 1)*length(conditions) + 3) = perm_cond(sequence(b), :)';
    end
    
    % save as .mat file
    filename = sprintf('%s_%s_stim_order.mat', study, subj_name);
    save(filename, 'stim_order')
    
    % append to global variable
    stim_order_all = cat(2, stim_order_all, stim_order);
end

% save the global cell
filename = sprintf('%s_all_stim_order.mat', study);
save(filename, 'stim_order_all')
clear perm_cond s subj_name b sequence stim_order filename

%% functions
function [fileID] = initializeMagProFilebi(filename,NofTrials)
% initializeMagProFile opens and writes in a text file specified by the
% filename variable, the first bloc of instructions required by the MagPro
% in order to read the experimental protocol file.
% NOTE : the extension given in the filename should be CG3!!!
fileID = fopen(filename,'w');
fprintf(fileID, '[Model Option]\r\n');
fprintf(fileID, 'ModelOption=3\r\n');
fprintf(fileID, '[Main Menu]\r\n');
fprintf(fileID, 'Mode=0\r\n');
fprintf(fileID, 'Current Direction=0\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Inter Pulse Interval=100\r\n');
fprintf(fileID, 'Burst Pulses=0\r\n');
fprintf(fileID, 'Pulse BA Ratio=100\r\n');
fprintf(fileID, '[Timing Menu]\r\n');
fprintf(fileID, 'Timing Control=0\r\n');
fprintf(fileID, 'Rep Rate=3\r\n');
fprintf(fileID, 'Pulses in train=1\r\n');
fprintf(fileID, 'Number of Trains=1\r\n');
fprintf(fileID, 'Inter Train Interval=30\r\n');
fprintf(fileID, '[Trigger Menu]\r\n');
fprintf(fileID, 'Trig Output=1\r\n');
fprintf(fileID, 'Twin Trig output=1\r\n');
fprintf(fileID, 'Twin Trig Input=0\r\n');
fprintf(fileID, 'Polarity Input=1\r\n');
fprintf(fileID, 'Polarity output=1\r\n');
fprintf(fileID, 'Delay Input Trig=0\r\n');
fprintf(fileID, 'Delay Output Trig=0\r\n');
fprintf(fileID, '[Configuration Menu]\r\n');
fprintf(fileID, 'Charge Delay=300\r\n');
fprintf(fileID, 'Auto Discharge Time=60\r\n');
fprintf(fileID, 'Prior Warning Sound=0\r\n');
fprintf(fileID, 'Coil Type Display=1\r\n');
fprintf(fileID, '[MEP Menu]\r\n');
fprintf(fileID, 'Time=5 ms/div\r\n');
fprintf(fileID, 'Sensitivity=500 uV/div\r\n');
fprintf(fileID, 'Scroll=0 ms\r\n');
fprintf(fileID, 'Curve No=0\r\n');
fprintf(fileID, 'Baseline=1\r\n');
fprintf(fileID, 'Lower Freq=100 Hz\r\n');
fprintf(fileID, 'Upper Freq=5 kHz\r\n');
fprintf(fileID, 'Trigger mode=0\r\n');
fprintf(fileID, 'Size=0\r\n');
fprintf(fileID, 'On Top=1\r\n');
fprintf(fileID, '[Protocol Setup]\r\n');
fprintf(fileID,'Number of Lines=%d\r\n',NofTrials);
fclose(fileID);
end
function [fileID] = writeMagPro_singlepulsebi(filename, trial_number, delay, amplitude)
fileID = fopen(filename,'a');
fprintf(fileID, '[protocol Line %d]\r\n',trial_number);
fprintf(fileID, 'Delay=%d\r\n',delay);
fprintf(fileID, 'Amplitude A Gain=%d\r\n',amplitude/10);
fprintf(fileID, 'Mode=0\r\n');
fprintf(fileID, 'Current Direction=1\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Burst Pulses=2\r\n');
fprintf(fileID, 'Inter Pulse Interval=10000\r\n');
fprintf(fileID, 'BA Ratio=100\r\n');
fprintf(fileID, 'Repetition Rate=10\r\n');
fprintf(fileID, 'Train Pulses=1\r\n');
fclose(fileID);
end
function [fileID] = MagPro_paired_bi_R(filename, trial_number, ITI, amplitude_CS, amplitude_TS, ISI)
fileID = fopen(filename,'a');
BAratio = round((amplitude_TS/amplitude_CS)*100);
if mod(BAratio,5) > 0
    if mod(BAratio,5) < 3
        BAratio = BAratio - mod(BAratio,5);
    end
    if mod(BAratio,5) >= 3
        BAratio = BAratio + (5 - mod(BAratio,5));
    end
end 
BAratio = round((amplitude_TS/amplitude_CS)*100);
fprintf(fileID, '[protocol Line %d]\r\n',trial_number);
fprintf(fileID, 'Delay=%d\r\n',ITI);
fprintf(fileID, 'Amplitude A Gain=%d\r\n',amplitude_CS/10);
fprintf(fileID, 'Mode=2\r\n');
fprintf(fileID, 'Current Direction=1\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Burst Pulses=2\r\n');
fprintf(fileID, 'Inter Pulse Interval=%d\r\n',ISI);
fprintf(fileID, 'BA Ratio=%d\r\n', BAratio);
fprintf(fileID, 'Repetition Rate=10\r\n');
fprintf(fileID, 'Train Pulses=1\r\n');
fclose(fileID)
end

                


