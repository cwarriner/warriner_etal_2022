function Process_Align_Dataset_2P_git(CNMF_filename,csv_filename,rabies_filename)

% Substantial portions written by Andrew Miri
clear;
plotfigs = 1;
saveFigs = 1;
processEMG = 1;

% Specify the onset & offset of EMG & encoder, if desired
% on = 136;
% off = 138; 
on = 156;
off = 158; 
% on = 38;
% off = 40;

%% Data Directories
    figDir = '<your_dir>';
    homeDir = 'C:\Users\user\Dropbox\MATLAB';
    signalsDir = 'C:\DATA\2P\CNMF Files';
    VoltageDir = 'C:\DATA\2P\Voltage Files';
    saveDir = 'C:\DATA\2P\CNMF Files';
    analyzDir = 'C:\DATA\2P\Analyzed Data';
    rabies_red_dir = 'C:\DATA\2P\Images_Red';
    mov_dir = 'C:\DATA\2P\Motion Corrected Data';

% %% Specify file names for bypass

% CNMF_filename = 'TSeries-04052018-Rd1_TR-001_CNMF_results';
% csv_filename = 'TSeries-04052018-Rd1_TR-001_Cycle00001_VoltageRecording_001.csv';
% rabies_img = 'SingleImage-04132018-Rd1_TR-006_Cycle00001_Ch2_000001.ome.tif';
% mov_filename = 'TSeries-04052018-Rd1_TR-001_mov_ch3_mov_nonrigid_corr';

CNMF_filename = 'TSeries-04092018-Rd1_TL-001_CNMF_results.mat';
csv_filename = 'TSeries-04092018-Rd1_TL-001_Cycle00001_VoltageRecording_001.csv';
rabies_img = 'SingleImage-04132018-Rd1_TR-006_Cycle00001_Ch2_000001.ome.tif';
mov_filename = 'TSeries-04092018-Rd1_TL-001_mov_ch3_mov_nonrigid_corr';

% CNMF_filename = 'TSeries-04092018-Rd1_TL-001_CNMF_results_p0.mat';
% csv_filename = 'TSeries-04092018-Rd1_TL-001_Cycle00001_VoltageRecording_001.csv';
% rabies_img = 'SingleImage-04132018-Rd1_TL-010_Cycle00001_Ch2_000001.ome.tif';
% mov_filename = 'TSeries-04092018-Rd1_TL-001_mov_ch3_mov_nonrigid_corr';
% 

%% File names
% if nargin < 3 || isempty(rabies_filename)
%     [rabies_filename, ~] = uigetfile('.tif','Select rabies image',rabies_red_dir);
% end
% 
% if nargin < 2 || isempty(voltage_loadname)
%     [voltage_loadname, ~] = uigetfile({'*.*'},'Select voltage file',VoltageDir);
% end
% 

    voltages_filename = strrep(CNMF_filename,'CNMF_results','voltages_ds'); % Rename CSV file for saving downsampled data

% 
% if nargin < 1 || isempty(CNMF_filename)
%     [CNMF_filename, ~] = uigetfile('*.mat','Select CNMF file',signalsDir);
% end

%% Parameters
frames_averaged = 4;
num_rows_orig = 256;
Y_galvo_col = 2; %column in csv of Y_galvo signal
Pockels_col = 3; %column in csv of Pockels blanking signal -IF NO POCKELS SAVED, COMMENT OUT
encoder_col = 4; %column in csv of encoder signal
solenoid_col = 5; %column in csv of solenoid signal
EMG_cols = 6:9; %columns of EMG signals in csv; change if one or more channels is dead
numIters_corr_sig = 1000;

samp_rate = 10000; %original sampling rate. Probably 10000
lp = 40; %low pass filter freq (Hz) for EMG filtering
smoFiltSTD = 10; %STD of Gaussian for EMG filtering, in ms
sub_rate = 10; %rate at which to downsample data - to go from 10K to 1k, use 10
Bi_ch = 1; %EMG header position where the biceps signal comes from
Tri_ch = 2; %EMG header position where the triceps signal comes from
EDC_ch = 3; %EMG header position where the EDC signal comes from
PL_ch = 4; %EMG header position where the PL signal comes from

alt_event_dur = 0.5; % duration of alternation event epochs, in sec
init_event_pre_dur = 0.3; %duration of quiesence before initiation events
init_event_post_dur = 1; %duration after initiation events
inact_dur = 0.4; %in s, duration of inactivity segments
modal_event_alt_size = 40;  % The X most similar alternation trials to keep
modal_event_init_size = 40;  % ? The X most similar initiation trials to keep
inact_event_size = 30;
var_thresh = 0.9;
merge_corr_thr = 0.8;
alpha_base = 0.05;

co1 = repmat([0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840],10,1);
%% Import CSV Data

    csv_readtime = tic;
    disp('LOADING CSV')
    Mat = csvread(fullfile(VoltageDir,csv_filename),1,0) ;
    % M = csvread(filename,R1,C1,[R1 C1 R2 C2]) 
    toc(csv_readtime)
    % 51 seconds for 20 minute movie

    %Downsample and separate signals
    if exist('Pockels_col','var')
        Pockels = downsample(Mat(:,Pockels_col),sub_rate);
    end
    Y_galvo_full = Mat(:,Y_galvo_col);
    Y_galvo = downsample(Mat(:,Y_galvo_col),sub_rate);
    encoder = downsample(Mat(:,encoder_col),sub_rate);
    solenoid = downsample(Mat(:,solenoid_col),sub_rate);
    EMG = downsample(Mat(:,EMG_cols),sub_rate);
    
    %save downsampled csv stuff for easier reloading
    if exist('Pockels_col','var')
        save(voltages_filename,'Pockels','Y_galvo','Y_galvo_full','encoder','solenoid','EMG')
    else
        save(voltages_filename,'Y_galvo','Y_galvo_full','encoder','solenoid','EMG')
    end

   % [~,peaks_y_galvo] = findpeaks(Y_galvo_full);
%     figure; hold on
%     plot(Y_galvo_full);
%     scatter(peaks_y_galvo,Y_galvo_full(peaks_y_galvo));
%     ylim([-10 10]); xlim([0 1000]);

%% Import CNMF results

load(fullfile(signalsDir,CNMF_filename));

num_cells = size(A2,2);
num_rows = num_rows_orig; %- length(rows_elim)
num_cols = num_rows_orig; %- length(cols_elim)
num_frames = size(C2,2);
num_muscles = size(EMG,2);

%Find the end of the last collected frame to set as the max time for event detection
a = Y_galvo > 0;
b = diff(a);
frame_samps = find(b==-1);

last_collected_frame_samp = frame_samps(num_frames*frames_averaged); % Error here: some misalignment or missing data..

c = unique(diff(frame_samps));
[d,ind] = min(Y_galvo(last_collected_frame_samp:last_collected_frame_samp+c(1)));
samp_end = last_collected_frame_samp+ind-1;
%pause

if exist('Pockels_col','var')
    Pockels = Pockels(1:samp_end);
end
Y_galvo = Y_galvo(1:samp_end);
encoder = encoder(1:samp_end);
solenoid = solenoid(1:samp_end);
EMG = EMG(1:samp_end,:);

%% Condition EMG
if processEMG
    EMG_smo = 0*EMG;
    for lag_idx = 1:size(EMG,2)
        EMG_ch_time = tic;
        EMG_smo(:,lag_idx) = filterEMG(EMG(:,lag_idx),lp,smoFiltSTD);
        toc(EMG_ch_time)
    end
    
    %Subtract the modal value off EMG signals to make baseline 0 - use last 5 minutes, when there is often not much activity
    for lag_idx = 1:size(EMG,2)
        num_datapoints = 100*samp_rate/sub_rate;
        [x,n] = hist(EMG_smo(end-num_datapoints:end,lag_idx),500);
        [m,ind] = max(x);
        signalMode = n(ind);
        EMG_smo(:,lag_idx) = EMG_smo(:,lag_idx) - signalMode;
    end
    
    %Normalize the EMG values to their 99.5th% value
    
    figure; hold on
    for lag_idx = 1:size(EMG,2)
        m = quantile(EMG_smo(:,lag_idx),0.995);
        EMG_smo(:,lag_idx) = EMG_smo(:,lag_idx)/m;
        %figure(next_fig_num); next_fig_num = next_fig_num+1; plot(EMG_smo(:,i))
        plot(EMG_smo(:,lag_idx),'LineWidth',2)
    end
end

%% Process Encoder
encoder_corrected = encoder(2:end,:);
diffEncoder = diff(encoder(:,1));
clockJumpIdx2 = find(diffEncoder < -1); % Find indices where the wheel made a complete clockwise turn & V jumped (-)
antiClockJumpIdx2 = find(diffEncoder > 1); % Find indices where the wheel made a complete anticlockwise turn & V jumped (+)
if ~isempty(clockJumpIdx2)
    for aa = 1:length(clockJumpIdx2)
        encoder_corrected(clockJumpIdx2(aa):end,1) = encoder_corrected(clockJumpIdx2(aa):end,1) + encoder(clockJumpIdx2(aa),1); % Add pre-jump V to rightward points
    end
end
if ~isempty(antiClockJumpIdx2)
    for bb = 1:length(antiClockJumpIdx2)
        encoder_corrected(antiClockJumpIdx2(bb):end,1) = encoder_corrected(antiClockJumpIdx2(bb):end,1) - ...
            ( encoder_corrected(antiClockJumpIdx2(bb),1) - encoder_corrected(antiClockJumpIdx2(bb)-1,1) ); % Subtract pre-jump V to rightward points
    end
end
% Filter Encoder w/ same parameters as EMG
    encoder_vel = diff(filterGauss(encoder_corrected,smoFiltSTD))';
    %figure; hold on;plot(encoder_vel); plot(encoder_vel_filt)


% %% Process Licks
% % Process ch 3 (otherwise Pockels) as licks if the voltage exceeds thresh
% if range(Pockels) > 3
%     licks_raw = Pockels;
%     licks_on = diff(licks_raw) > 4; % remove brief V jumps
% else
%     no_lick_info = 1;
% end

%% Process Reward Solenoid
% Find onset time of reward trains and mark the pulse train onset index with the # of pulses in the train

rew_on = diff(solenoid) > 4;                                                % Find reward onsets
rew_info = zeros(length(rew_on),1);                                         % Set up vector for reporting pulse train identification

samp_rate_dwnsmp = samp_rate/sub_rate;
rew_time = linspace(0,length(rew_on)/samp_rate_dwnsmp,length(rew_on));      % Make time
fwd_window = 0.5*samp_rate_dwnsmp;                                          % Look ahead for solenoid clicks that belong to a reward train ** up to FR5 compatible

ss = 1;
while ss <= length(rew_on) - fwd_window
    pulse_length = sum(rew_on(ss:ss + fwd_window)); 
    if rew_on(ss) ~= 0                                                      % Only mark pulse trains that share an onset w/the current index
        rew_info(ss) = pulse_length;                                        % Mark the length of the pulse train
        if pulse_length >= 1
            ss = ss + fwd_window;                                           % Skip ahead by forward looking window value
        else
            ss = ss + 1;                                                    % Or move on to next index
        end
    else
        ss = ss + 1;                                                        % Or move on to next index
    end
end
    
% figure; hold on
% plot(rew_time(1:end-1),encoder_vel')
% plot(rew_time,rew_on/100)
% plot(rew_time,rew_info)
% disp('Check for accuracy')
% pause

%% Find peaks in EMG

pre_win_alt = floor(alt_event_dur*samp_rate/sub_rate/2);

if processEMG
    % Peaks in Triceps
    [~, locs_tri] = findpeaks(EMG_smo(:,Tri_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_tri(locs_tri < pre_win_alt) = [];
    locs_tri(locs_tri > size(EMG,1) - pre_win_alt) = [];
    
    % Peaks in Biceps
    [~, locs_bi] = findpeaks(EMG_smo(:,Bi_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_bi(locs_bi < pre_win_alt) = [];
    locs_bi(locs_bi > size(EMG,1) - pre_win_alt) = [];
    
    % Peaks in EDC
    [~, locs_EDC] = findpeaks(EMG_smo(:,EDC_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_EDC(locs_EDC < pre_win_alt) = [];
    locs_EDC(locs_EDC > size(EMG,1) - pre_win_alt) = [];
    
    % Peaks in PL
    [~, locs_PL] = findpeaks(EMG_smo(:,PL_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_PL(locs_PL < pre_win_alt) = [];
    locs_PL(locs_PL > size(EMG,1) - pre_win_alt) = [];
    
    
    %% Find alternation events (Bi/Tri)
    % Find modal duration from tri peak to next bi peak
    durs_pos = zeros(length(locs_tri),1);
    durs_neg = zeros(length(locs_tri),1);
    for lag_idx = 1:length(locs_tri) % For each tri peak index...
        diffs = (locs_bi - locs_tri(lag_idx)); % find the time shifts between all bi peaks and that tri peak
        diffs_pos = diffs(diffs>0); % separate the bi peaks that come after or before the tri peak
        diffs_neg = diffs(diffs<0);
        if ~isempty(diffs_pos)
            durs_pos(lag_idx) = min(diffs_pos); % and find the first bi peak that comes after the tri peak
        end
        if ~isempty(diffs_neg)
            durs_neg(lag_idx) = max(diffs_neg); % and find the last bi peak that comes before the tri peak
        end
    end
    
    [n,bins_pos] = histc(durs_pos,50:5:pre_win_alt); % Bin the post-tri bi peaks; peaks must be at least 50 ms later
    [m,ind_pos] = max(n);   % The peak in histogram corresponds to the mode
    dur_pos_mode = 50+ind_pos*5-2.5;
    median_dur_pos = median(durs_pos);
    
    [n,bins_neg] = histc(durs_neg,-pre_win_alt:5:-50); % Bin the pre-tri bi peaks
    [m,ind_neg] = max(n);
    dur_neg_mode = -pre_win_alt+ind_neg*5-2.5;
    median_dur_neg = median(durs_neg);
    
    modal_events_alt = [];
    width = 0;
    while length(modal_events_alt) < modal_event_alt_size   % # of alternation events determined by modal_event_alt_size
        modal_pos_events = find(bins_pos <= ind_pos+width & bins_pos >= ind_pos-width);     % Find post-tri bi peaks that are within the window of mode+width and mode-width
        modal_neg_events = find(bins_neg <= ind_neg+width & bins_neg >= ind_neg-width);     % Repeat for pre-tri bi peaks
        modal_events_alt = intersect(modal_pos_events,modal_neg_events);    % Find the indices where they overlap
        width = width+1;    % Expand the width of the window... until enough alternation events have been found
    end
    
    num_alt_events = length(modal_events_alt);
    events_alt_emg = zeros(num_alt_events,pre_win_alt*2+1,4);
    for lag_idx = 1:num_alt_events % Get EMG traces that correspond to each alternation event
        events_alt_emg(lag_idx,:,1) = EMG_smo(locs_tri(modal_events_alt(lag_idx))-pre_win_alt:locs_tri(modal_events_alt(lag_idx))+pre_win_alt,Bi_ch); % For bi
        events_alt_emg(lag_idx,:,2) = EMG_smo(locs_tri(modal_events_alt(lag_idx))-pre_win_alt:locs_tri(modal_events_alt(lag_idx))+pre_win_alt,Tri_ch); % Tri
        events_alt_emg(lag_idx,:,3) = EMG_smo(locs_tri(modal_events_alt(lag_idx))-pre_win_alt:locs_tri(modal_events_alt(lag_idx))+pre_win_alt,EDC_ch); % EDC
        events_alt_emg(lag_idx,:,4) = EMG_smo(locs_tri(modal_events_alt(lag_idx))-pre_win_alt:locs_tri(modal_events_alt(lag_idx))+pre_win_alt,PL_ch); % PL
        
    end
    
    % Plot modal alt averages
    figure;
    time_alt = -pre_win_alt:pre_win_alt;
    boundedline(time_alt,mean(events_alt_emg(:,:,2)),std(events_alt_emg(:,:,2))/sqrt(num_alt_events),'-b','transparency', 0.5)
    hold on;
    boundedline(time_alt,mean(events_alt_emg(:,:,1)),std(events_alt_emg(:,:,1))/sqrt(num_alt_events),'-r','transparency', 0.5)
    
    xlabel('Time from triceps peak (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Alternation EMG')
    ylim([-0.05 1])
    hold off
    
    %% Find EMG initiation events
    
    time = 0.001:0.001:0.001*size(EMG_smo,1);
    pre_win_init = init_event_pre_dur*samp_rate/sub_rate;
    post_win_init = init_event_post_dur*samp_rate/sub_rate;
    figure;plot(time,EMG_smo(:,Bi_ch))  % plot EMG to find baseline
    if isempty(on)
        bi_baseline_segment_begin = input('Baseline bi segment begins? in sec  ');     % Choose time periods for determining baseline
        bi_baseline_segment_end = input('Baseline bi segment ends? in sec  ');
    else
        bi_baseline_segment_begin = on;
        bi_baseline_segment_end = off;
    end
    bi_baseline_std = std(EMG_smo(bi_baseline_segment_begin*samp_rate/sub_rate:bi_baseline_segment_end*samp_rate/sub_rate,Bi_ch));  % Find STD of baseline period
    bi_activation_th = 10*bi_baseline_std;   % Threshold for crossings: 10x STD of baseline
    dummy = diff(EMG_smo(:,Bi_ch) > bi_activation_th);   % Find crossings
    crossings = find(dummy == 1);
    crossings(crossings<pre_win_init) = [];   % Get ride of out-of-bounds crossings
    crossings(crossings>length(EMG_smo(:,Bi_ch))-post_win_init) = [];  % Get ride of out-of-bounds crossings
    figure; plot(time,EMG_smo(:,Tri_ch))  % Repeat for Tri
    if isempty(on)
        tri_baseline_segment_begin = input('Baseline tri segment begins? in sec  ');
        tri_baseline_segment_end = input('Baseline tri segment ends? in sec  ');
    else
        tri_baseline_segment_begin = on;
        tri_baseline_segment_end = off;
    end
    tri_baseline_std = std(EMG_smo(tri_baseline_segment_begin*samp_rate/sub_rate:tri_baseline_segment_end*samp_rate/sub_rate,Tri_ch));
    tri_activation_th = 10*tri_baseline_std;
    
    elim = 1;
    for lag_idx = 2:length(crossings)
        if sum(EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)-1,Bi_ch) > bi_activation_th) > 0 % If there is persistent activity w/in 300msec of a crossing, eliminate
            elim = [elim; lag_idx];
        elseif sum(EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)-1,Tri_ch) > tri_activation_th) > 0 % In either channel (Tri, here)
            elim = [elim; lag_idx];
        end
    end
    crossings(elim) = [];
    
    %Find modal init event with peaks
    durs_bi = zeros(length(crossings),1);
    durs_tri = zeros(length(crossings),1);
    durs_EDC = zeros(length(crossings),1);
    durs_PL = zeros(length(crossings),1);
    
    for lag_idx = 1:length(crossings)
        diffs = (locs_bi - crossings(lag_idx));   % Find the distance from each bi peak to each bi threshold crossing
        diffs_bi = diffs(diffs>0);	% Peak must come after the crossing
        if ~isempty(diffs_bi)
            durs_bi(lag_idx) = min(diffs_bi);	% Get the first peak to come after the threshold crossing; result is the durations between threshold crossings and bi peaks
        end
        diffs = (locs_tri - crossings(lag_idx));	% Repeat for tri
        diffs_tri = diffs(diffs>0);
        if ~isempty(diffs_tri)
            durs_tri(lag_idx) = min(diffs_tri);   % The result is the durations between threshold crossings and tri peaks
        end
        diffs = (locs_EDC - crossings(lag_idx));	% Repeat for EDC
        diffs_EDC = diffs(diffs>0);
        if ~isempty(diffs_EDC)
            durs_tri(lag_idx) = min(diffs_EDC);   % Time from crossings to EDC peak
        end
        diffs = (locs_PL - crossings(lag_idx));	% Repeat for PL
        diffs_PL = diffs(diffs>0);
        if ~isempty(diffs_PL)
            durs_PL(lag_idx) = min(diffs_PL);     % Time from crossings to PL peak
        end
    end
    
    [n,bins_bi] = histc(durs_bi,0:5:post_win_init/2);   % Bin the vector of times from crossings to peaks, bi
    [m,ind_bi] = max(n);
    dur_bi_mode = ind_bi*5-2.5;
    
    [n,bins_tri] = histc(durs_tri,2*dur_bi_mode:5:post_win_init); %peaks must be after bi peaks : ??????
    [m,ind_tri] = max(n);
    dur_tri_mode = 2*dur_bi_mode+ind_tri*5-2.5;
    
    modal_events_init = [];
    width = 0;
    while length(modal_events_init) < 2 % why 2?
        modal_bi_events = find(bins_bi <= ind_bi+width & bins_bi >= ind_bi-width);
        modal_tri_events = find(bins_tri <= ind_tri+width & bins_tri >= ind_tri-width);
        modal_events_init = intersect(modal_bi_events,modal_tri_events);
        width = width+1;
    end
    
    %Now find the modal time series and compute distance from it.
    crossings_ts = zeros(length(crossings),pre_win_init+post_win_init+1,4);
    for lag_idx = 1:length(crossings)
        crossings_ts(lag_idx,:,1) = EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,Bi_ch)';
        crossings_ts(lag_idx,:,2) = EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,Tri_ch)';
        crossings_ts(lag_idx,:,3) = EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,EDC_ch)';
        crossings_ts(lag_idx,:,4) = EMG_smo(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,PL_ch)';
    end
    
    crossings_ts_clique = [ [mean(crossings_ts(modal_events_init,pre_win_init+1:end,1)) mean(crossings_ts(modal_events_init,pre_win_init+1:end,2))];...
        [crossings_ts(:,pre_win_init+1:end,1) crossings_ts(:,pre_win_init+1:end,2)] ]; % Looks like top row = mean of modal initiation events for bi and tri, concatenated.  All other rows = all threshold crossings for bi/tri, concatenated
    dists = squareform(pdist(crossings_ts_clique)); % Find pairwise distand between each crossings instance: probably used to find the most similar trials
    dists = dists(2:end,1); %the top row gives the distance of each row from the top row of the input matrix, which had the modal event mean as the top row.
    [sorted,sorted_inds] = sort(dists);
    if length(sorted_inds) < modal_event_init_size
        modal_event_init_size = length(sorted_inds);
    end
    keep_init = sorted_inds(1:modal_event_init_size); % Keep X most similar trials?
    
    events_init_emg(:,:,1) = crossings_ts(keep_init,:,1);
    events_init_emg(:,:,2) = crossings_ts(keep_init,:,2);
    events_init_emg(:,:,3) = crossings_ts(keep_init,:,3);
    events_init_emg(:,:,4) = crossings_ts(keep_init,:,4);
    num_init_events = length(keep_init);
    
    %Plot modal init averages
    figure;
    time_init = -pre_win_init:post_win_init;
    boundedline(time_init,mean(events_init_emg(:,:,2)),std(events_init_emg(:,:,2))/sqrt(num_init_events),'-b','transparency', 0.5,'alpha')
    hold on;
    boundedline(time_init,mean(events_init_emg(:,:,1)),std(events_init_emg(:,:,1))/sqrt(num_init_events),'-r','transparency', 0.5,'alpha')
    boundedline(time_init,mean(events_init_emg(:,:,3)),std(events_init_emg(:,:,3))/sqrt(num_init_events),'-m','transparency', 0.5,'alpha')
    boundedline(time_init,mean(events_init_emg(:,:,4)),std(events_init_emg(:,:,4))/sqrt(num_init_events),'-c','transparency', 0.5,'alpha')
    xlabel('Time from biceps threshold crossing (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Initiation EMG')
    ylim([-0.05 1])
    hold off
    
    % Find inactivity periods
%     inact_dur_samps = inact_dur*samp_rate/sub_rate;
%     summed_EMG = [];
%     wb=waitbar(0,'Finding inactive segments');
%     for t = 1:10:size(EMG_smo,1)-inact_dur_samps % Sum the EMG segment by segment
%         summed_EMG = [summed_EMG; t sum(sum(EMG_smo(t:t+inact_dur_samps-1,:)))];
%         h = waitbar(t/size(EMG_smo,1),wb,'Finding inactive segments');
%     end
%     close(h)
%     summed_EMG_sorted = sortrows(summed_EMG,2); % Sort smallest to largest
%     
%     %Now eliminate the overlapping chunks..????????????????
%     viols = 1;
%     while ~isempty(viols)
%         dummy = abs(bsxfun(@minus,summed_EMG_sorted(1:inact_event_size,1),summed_EMG_sorted(1:inact_event_size,1)')) < inact_dur_samps; % wut
%         [row,col] = find(triu(dummy,1));
%         viols = unique(col);
%         summed_EMG_sorted(viols,:) = [];
%     end
%     
%     events_inact_emg = zeros(inact_event_size,inact_dur_samps,4);
%     for lag_idx = 1:inact_event_size
%         events_inact_emg(lag_idx,:,1) = EMG_smo(summed_EMG_sorted(lag_idx,1):summed_EMG_sorted(lag_idx,1)+inact_dur_samps-1,Bi_ch);
%         events_inact_emg(lag_idx,:,2) = EMG_smo(summed_EMG_sorted(lag_idx,1):summed_EMG_sorted(lag_idx,1)+inact_dur_samps-1,Tri_ch);
%         events_inact_emg(lag_idx,:,3) = EMG_smo(summed_EMG_sorted(lag_idx,1):summed_EMG_sorted(lag_idx,1)+inact_dur_samps-1,EDC_ch);
%         events_inact_emg(lag_idx,:,4) = EMG_smo(summed_EMG_sorted(lag_idx,1):summed_EMG_sorted(lag_idx,1)+inact_dur_samps-1,PL_ch);
%     end
%     num_inact_events = inact_event_size;
%     
    % %Plot inact averages
    % figure;
    % time_inact = 1:inact_dur_samps;
    % boundedline(time_inact,mean(events_inact_emg(:,:,2)),std(events_inact_emg(:,:,2))/sqrt(num_inact_events),'-b','transparency', 0.5)
    % hold on;
    % boundedline(time_inact,mean(events_inact_emg(:,:,1)),std(events_inact_emg(:,:,1))/sqrt(num_inact_events),'-r','transparency', 0.5)
    % boundedline(time_inact,mean(events_inact_emg(:,:,3)),std(events_inact_emg(:,:,3))/sqrt(num_inact_events),'-m','transparency', 0.5)
    % boundedline(time_inact,mean(events_inact_emg(:,:,4)),std(events_inact_emg(:,:,4))/sqrt(num_inact_events),'-c','transparency', 0.5)
    % xlabel('Time (ms)')
    % ylabel('Normalized EMG magnitude (ms)')
    % title('Inactivation EMG')
    % ylim([-0.05 1])
    % hold off
end
%% Process and align encoder 

% Peaks in Encoder
encoder_vel_norm = (encoder_vel - min(encoder_vel)) ./ (max(encoder_vel) - min(encoder_vel));
[~, locs_enc] = findpeaks(encoder_vel_norm,'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
locs_enc(locs_enc < pre_win_alt) = [];
locs_enc(locs_enc > size(EMG,1) - pre_win_alt) = [];
time = 0.001:0.001:0.001*size(encoder_vel_norm,1);
figure; plot(time,encoder_vel_norm); hold on; scatter(time(locs_enc),encoder_vel_norm(locs_enc));
set(gca,'FontSize',12)
title('Encoder Velocity & Detected Peaks','FontSize',14);
xlabel('Time (sec)','FontSize',14)
ylabel('Normalized Velocity','FontSize',14)


pre_win_init = init_event_pre_dur*samp_rate/sub_rate;
post_win_init = init_event_post_dur*samp_rate/sub_rate;
if isempty(on)
    end_baseline_segment_begin = input('Baseline segment begins? in s  ');
    enc_baseline_segment_end = input('Baseline segment ends? in s  ');
else
    end_baseline_segment_begin = on;
    enc_baseline_segment_end = off;
end
encoder_vel_norm_sub = encoder_vel_norm - mean(encoder_vel_norm(end_baseline_segment_begin*samp_rate/sub_rate:enc_baseline_segment_end*samp_rate/sub_rate));
%encoder_vel_norm_sub = encoder_vel;
enc_baseline_std = std(encoder_vel_norm_sub(end_baseline_segment_begin*samp_rate/sub_rate:enc_baseline_segment_end*samp_rate/sub_rate));
enc_activation_th = 20*enc_baseline_std;
dummy = diff(encoder_vel_norm_sub > enc_activation_th);
crossings2 = find(dummy == 1);
crossings2(crossings2<pre_win_init) = [];
crossings2(crossings2>length(encoder_vel_norm_sub)-post_win_init) = [];
disp('Processing Encoder')

all_crossings2 = crossings2;
elim = 1;
for lag_idx = 2:length(crossings2) %get rid of those with other crossings (activity) within 300 ms
    if sum(encoder_vel_norm_sub(crossings2(lag_idx)-pre_win_init:crossings2(lag_idx)-1) > enc_activation_th) > 0
        elim = [elim; lag_idx];
    end
end
crossings2(elim) = [];

% Repurpose EMG code for encoder; Find modal init event with peaks
durs_enc = zeros(length(crossings2),1);
for lag_idx = 1:length(crossings2)
    diffs = (locs_enc - crossings2(lag_idx));
    diffs_enc = diffs(diffs>0);
    if ~isempty(diffs_enc)
        durs_enc(lag_idx) = min(diffs_enc);
    end
end
[n,bins_enc] = histc(durs_enc,0:5:post_win_init/2); 
[m,ind_enc] = max(n);
dur_enc_mode = ind_enc*5-2.5;

%Now find the modal time series and compute distance from it.
crossings_enc = zeros(length(crossings2),pre_win_init+post_win_init+1);
for lag_idx = 1:length(crossings2)
    crossings_enc(lag_idx,:) = encoder_vel_norm_sub(crossings2(lag_idx)-pre_win_init:crossings2(lag_idx)+post_win_init)';
end

figure;
time_init = -pre_win_init:post_win_init;
boundedline(time_init,mean(crossings_enc(:,:)),std(crossings_enc(:,:))/sqrt(size(crossings_enc,1)),'-b','transparency', 0.5,'alpha')
set(gca,'FontSize',12)
xlabel('Time from encoder threshold crossing (ms)','FontSize',14)
ylabel('Normalized encoder magnitude','FontSize',14)
title('Initiation Encoder','FontSize',14)
if saveFigs
    cd(figDir)
    print -painters -depsc Encoder_sequence_binned_durs.eps 

end
hold off
figure; hold on
for oo = 1:length(crossings2)
    plot(time_init,crossings_enc(oo,:),'color',[.3 .3 .3])
end

%% Condition Ca data
% 1) Identify the temporal bounds of Ca aquisition 
% 2) Gaussian filter Ca data with same parameters as EMG
% 3) Remove ROIs that are near the periphery of the FOV
% 4) Average highly correlated Ca signals & remove all but one of the constituent ROIs
% 5) Align Ca activity to the initiation of biceps activity
% 6) Align Ca activity to initiation of encoder activity
% 7) Sort out bouts of encoder activity of different lengths
% 8) Align Ca acitivity to different bout durations of encoder activity

% All taken from Miri's script
disp('Coniditioning Ca Data')
galvo_bnds = quantile(Y_galvo_full(end/4:3*end/4),[0 1]); %avoid the ends where Y_galvo may hit extremes of range not seen during scanning at higher zone
Cell_ts_raw = zeros(size(Y_galvo,1),num_cells);
Cell_ts_ca = zeros(size(Y_galvo,1),num_cells);
Cell_ts_sp = zeros(size(Y_galvo,1),num_cells);
%raw_rois = raw_rois';
for lag_idx = 1:num_cells
    %Mat2 = reshape(A_or(:,lag_idx),num_rows,num_cols,1); % Which one?  I
    % C_df is now from unordered componenets
    Mat2 = reshape(A2(:,lag_idx),num_rows,num_cols,1); 

    cOM = centerOfMass(Mat2);
    cell_row_pos = cOM(1);
    row_galvo_pos = galvo_bnds(2)- abs(diff(galvo_bnds)/num_rows)*cell_row_pos;
    a = Y_galvo_full > row_galvo_pos;
    b = diff(a);
    row_times_full = find(b==-1);
    row_times = round(row_times_full/sub_rate); %Get back to downsampled time points
    row_times = row_times(1:num_frames*frames_averaged);
    c = diff(row_times);
    d = unique(c);
    half_frame_time = floor(d(1)/2);
    for j = 1:num_frames
        bnds = [max([1 row_times(frames_averaged*(j-1)+1)-half_frame_time])  min([size(Y_galvo,1) row_times(frames_averaged*j)+half_frame_time])];
        %Cell_ts_raw(bnds(1):bnds(2),lag_idx) = raw_rois(lag_idx,j); % Don't
        %use? maybe raw_rois
        Cell_ts_ca(bnds(1):bnds(2),lag_idx) = C_df(lag_idx,j);
        %Cell_ts_ca(bnds(1):bnds(2),lag_idx) = C2(lag_idx,j);
        Cell_ts_sp(bnds(1):bnds(2),lag_idx) = S2(lag_idx,j)*frames_averaged; %I think this is right - sp is computed as spikes per time bin/frame?
    end
    
    %Smooth Cell time series to be equivalent to EMG, getting rid of rough edges
    %Cell_ts_raw(:,lag_idx) = filterGauss(Cell_ts_raw(:,lag_idx),smoFiltSTD);
    Cell_ts_ca(:,lag_idx) = filterGauss(Cell_ts_ca(:,lag_idx),smoFiltSTD);
    Cell_ts_sp(:,lag_idx) = filterGauss(Cell_ts_sp(:,lag_idx),smoFiltSTD);
    
end
%% Look at raw data

% figure;hold on
% for oo = 1:size(Cell_ts_ca,2)
%     plot(Cell_ts_ca(:,oo))
% end

%% Remove ROIs that are close to the periphery and may shift out of view during movement
% Use spatial footprints to exclude ROIs that fall at the edge of the FOV
figure;
%[Coor,json_file] = plot_contours(A2,Cn,options,1);
[Coor,~] = plot_contours(A_or,Cn,options,1);

border_exclude = 10; % ROIs with any pixels this close to the border are excluded
exclude_ROI_border = [];
for a = 1:length(Coor)
    if any(any(Coor{a} < border_exclude)) || any(any(Coor{a} > 256-border_exclude))
        exclude_ROI_border = [exclude_ROI_border a];
    end
end
%% Choose data type to analyze
real_Cell_ts_ca = Cell_ts_ca;
%real_Cell_ts_ca = Cell_ts_raw;
real_Cell_ts_ca(:,exclude_ROI_border) = nan;

%% Find highly similar ROIs using pairwise correlation coefficient
rho = triu(corr(real_Cell_ts_ca),1);
hi_corr = rho > 0.8; % Correlation greater than 0.8 will signal same cell
figure;imagesc(hi_corr)

pruned_Cell_ts_ca = NaN(length(Cell_ts_ca),num_cells);
pruned_Cell_ts_sp = Cell_ts_sp;
redunant_ROI=[];
for a = 1:length(hi_corr)
    if any(hi_corr(a,:))
        pruned_Cell_ts_ca(:,a) = mean([real_Cell_ts_ca(:,find(hi_corr(a,:))) real_Cell_ts_ca(:,a)], 2); % Replace current ROI activity with mean activity of all highly similar ROIs + activity from current ROI in loop
        redunant_ROI = [redunant_ROI find(hi_corr(a,:))]; % So these are the leftover highly correlated ROIs that should not contribute to the analysis
    elseif sum(~isnan(real_Cell_ts_ca(:,a))) ~= 0
        pruned_Cell_ts_ca(:,a) = real_Cell_ts_ca(:,a); % If there are no other ROIs with highly correlated activity (and it was not classified as a border ROI), just keep it
    end
end
redunant_ROI=unique(redunant_ROI);
pruned_Cell_ts_ca(:,redunant_ROI) = nan; % The highly correlated ROIs that went into an average are now made to be nan

ROI_keeper_IDs = find(~isnan(pruned_Cell_ts_ca(1,:))); % ID the ROIs that are valid

totrash = isnan(pruned_Cell_ts_ca(1,:));

pruned_Cell_ts_ca(:,totrash) = []; % And remove invalid ROIs for simplicity
pruned_Cell_ts_sp(:,totrash) = [];

%A2_pruned = A2(:,ROI_keeper_IDs); % Spatial footprints of the remaining ROIs
A2_pruned = A_or(:,ROI_keeper_IDs); % Spatial footprints of the remaining ROIs

figure;
[Coor_pruned] = plot_contours(A2_pruned,Cn,options,1);

    
%% Align Ca data to EMG
bi_crossings_Ca = zeros(length(crossings),pre_win_init+post_win_init+1,size(pruned_Cell_ts_ca,2));
bi_crossings_Sp = zeros(length(crossings),pre_win_init+post_win_init+1,size(pruned_Cell_ts_ca,2));

for lag_idx = 1:length(crossings)
    for cc = 1:size(pruned_Cell_ts_ca,2)
        bi_crossings_Ca(lag_idx,:,cc) = pruned_Cell_ts_ca(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,cc)';
        bi_crossings_Sp(lag_idx,:,cc) = pruned_Cell_ts_ca(crossings(lag_idx)-pre_win_init:crossings(lag_idx)+post_win_init,cc)';
    end
end
time_init = -pre_win_init:post_win_init;

% Z Score Ca activity
[Z_bi_Ca, mu_bi_Ca, sig_bi_Ca] = zscore(bi_crossings_Ca(:,:,:),1,1);
Z_avg_Ca = mean(Z_bi_Ca);
Z_sig_Ca = std(Z_bi_Ca);

% Z Score Sp activity
[Z_bi_Sp, mu_bi_Sp, sig_bi_Sp] = zscore(bi_crossings_Sp(:,:,:),1,1);
Z_avg_Sp = mean(Z_bi_Sp);
Z_sig_Sp = std(Z_bi_Sp);


% Plot all ROIs, different ways   
% num_init_events is +1 ????
numSubAxis = ceil(sqrt(size(pruned_Cell_ts_ca,2)));
if plotfigs
    figure;
    for cc = 1:size(pruned_Cell_ts_ca,2)
        h1(cc) = subplot(numSubAxis,numSubAxis,cc); hold on
        %boundedline(time_init,mean(bi_crossings_Ca(:,:,cc)),std(bi_crossings_Ca(:,:,cc))/sqrt(num_init_events),'-b','transparency', 0.1,'alpha')
        boundedline(time_init,mean(bi_crossings_Ca(:,:,cc))/max(mean(bi_crossings_Ca(:,:,cc))), std(bi_crossings_Ca(:,:,cc))/sqrt(num_init_events-1)/max(mean(bi_crossings_Ca(:,:,cc))),'-b','transparency', 0.1,'alpha')
        %boundedline(time_init,mean(bi_crossings_Sp(:,:,cc))/max(mean(bi_crossings_Sp(:,:,cc))), std(bi_crossings_Sp(:,:,cc))/sqrt(num_init_events)/max(mean(bi_crossings_Sp(:,:,cc))),'-c','transparency', 0.1,'alpha')
        
        %     for ii = 1:size(Cell_ts_ca,2)
        %        plot(time_init,Z_bi_Ca(ii,:,cc),'Color',[.3 .3 .3])
        %     end
        %boundedline(time_init,mu_bi_Ca(:,:,cc),sig_bi_Ca(:,:,cc)/sqrt(num_init_events),'-b','transparency', 0.05,'alpha')
        %boundedline(time_init,mu_bi_Sp(:,:,cc),Z_sig_Sp(:,:,cc)/sqrt(num_init_events),'-g','transparency', 0.05,'alpha')
        
        %ylabel('Normalized Ca activity','fontsize',9)
        ylabel('Ca activity [Z Score]','fontsize',9)
        xlabel('Time relative to onset (msec)','fontsize',9)
    end
    linkaxes(h1,'x')
    suptitle('Ca data aligned to bi EMG onset')
end
% figure; hold on
% for cc = 1:size(Cell_ts_ca,2)
%     %plot(time_init, mean_crossings_Ca_norm(:,:,cc));
%     %boundedline(time_init,mean(crossings_Ca(:,:,cc)),std(crossings_Ca(:,:,cc))/sqrt(num_init_events),'-b','transparency', 0.1,'alpha')
%     boundedline(time_init,mu_bi_Ca(:,:,cc),sig_bi_Ca(:,:,cc)/sqrt(num_init_events),'-b','transparency', 0.03,'alpha')
% end
% title('Ca data aligned to bi EMG onset')
% ylabel('Ca activity [Z Score]')
% xlabel('Time relative to onset (msec)')

%% Plot EMG & Encoder together
    figure; 
    h(1) = subtightplot(2,6,[1:6]); hold on;
    for ee = 1:size(EMG,2)
        m = quantile(EMG_smo(:,ee),0.995);
        EMG_smo(:,ee) = EMG_smo(:,ee)/m;
        %figure(next_fig_num); next_fig_num = next_fig_num+1; plot(EMG_smo(:,i))
        plot(linspace(0,length(EMG_smo(:,ee))/1000,length(EMG_smo(:,ee))),EMG_smo(:,ee),'LineWidth',1)
    end
        set(h(1),'xtick',[]);
        axes('Color','none','XColor','none');
    h(2) = subtightplot(2,6,[7:12])
    plot(time, encoder_vel)
    linkaxes(h,'x');
%% Align Ca data to encoder data
disp('Aligning Ca data to encoder onset')
enc_crossings_Ca = zeros(length(crossings2),pre_win_init+post_win_init+1,size(pruned_Cell_ts_ca,2));
for lag_idx = 1:length(crossings2)
    for cc = 1:size(pruned_Cell_ts_ca,2)
        enc_crossings_Ca(lag_idx,:,cc) = pruned_Cell_ts_ca(crossings2(lag_idx)-pre_win_init:crossings2(lag_idx)+post_win_init,cc)';
    end
end

% Plot Ca activity relative to all encoder activity (also EMG for comparison)
disp('Plotting Ca activity relative to EMG and encoder')
numSubAxis = ceil(sqrt(size(pruned_Cell_ts_ca,2)));
if plotfigs
    figure;
    for cc = 1:size(enc_crossings_Ca,3)
        h2(cc) = subplot(numSubAxis,numSubAxis,cc); hold on
        boundedline(time_init,mean(bi_crossings_Ca(:,:,cc))/max(mean(bi_crossings_Ca(:,:,cc))), std(bi_crossings_Ca(:,:,cc))/sqrt(num_init_events-1)/max(mean(bi_crossings_Ca(:,:,cc))),'-b','transparency', 0.1,'alpha')
        boundedline(time_init,mean(enc_crossings_Ca(:,:,cc))/max(mean(enc_crossings_Ca(:,:,cc))),std(enc_crossings_Ca(:,:,cc)/max(mean(enc_crossings_Ca(:,:,cc))))/sqrt(size(enc_crossings_Ca,3)-1),'-g','transparency', 0.1,'alpha')
        %boundedline(time_init,mean(enc_crossings_Ca(:,:,cc)),std(enc_crossings_Ca(:,:,cc))/sqrt(num_init_events),'-b','transparency', 0.1,'alpha')
    end
    linkaxes(h2,'x')
end
% Sort out different lengths of encoder activity
disp('Sorting encoder activity by duration')
fwd_window2 = 300;
enc_thresh = [];
enc_thresh = encoder_vel_norm_sub > enc_activation_th;                      % Find where encoder is higher than threshold
ss = 1;
while ss <= length(encoder_vel) - fwd_window2 
    if enc_thresh(ss) == 1 && sum(enc_thresh(ss:ss + fwd_window2)) > 1      % First combine encoder epochs that are within 'fwd_window2' of each other
        if isempty(intersect(ss+1, crossings2))                             % And don't overlap with detected crossings
            enc_thresh(ss + 1) = 1;
        end
    end
    ss = ss + 1;
end
if plotfigs
    figure; hold on % Check for accuracy
    plot(time,encoder_vel_norm_sub);
    scatter(time(crossings2),encoder_vel_norm_sub(crossings2))
    plot(time,enc_thresh,'k','LineWidth',1.5)
    set(gca,'FontSize',12)
    xlabel('Time (sec)','FontSize',14)
    ylabel('Normalized encoder magnitude','FontSize',14)
    title('Encoder sequence sorted','FontSize',14)
    if saveFigs
        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'Encoder_sequence_chunking_examp.eps'))
    end

end
% For each encoder timestamp measure pulse duration; Bin & Histogram
pulse_onsets = [];
pulse_offsets=[];
pulse_onsets = find(diff(enc_thresh) == 1);                                 % Find onsets
pulse_offsets = find(diff(enc_thresh) == -1);                               % Find offsets

if pulse_offsets(1) < pulse_onsets(1)                                       % Trim pulse offsets that don't have an onset (onset before aquisition starts)
    pulse_offsets(1) = [];
end

if pulse_onsets(end) > pulse_offsets(end)                                   % Trim threshold crossing events w/o an offset (offset after aquisition ends)
    pulse_onsets(end) = [];
end

crossings2_trim = crossings2;
pulse_onsets_trim = [];
pulse_offsets_trim = [];
if length(crossings2_trim) ~= length(pulse_offsets)
    for w = 1:length(pulse_onsets)
        if ~any(pulse_onsets(w) == crossings2)                              % If a crossing wasn't detected but there was a pulse
            pulse_onsets_trim = [pulse_onsets_trim; w];
            pulse_offsets_trim = [pulse_offsets_trim; w];
        end
    end
end
pulse_onsets(pulse_onsets_trim) = [];
pulse_offsets(pulse_offsets_trim) = [];

enc_pulse_lengths = pulse_offsets - crossings2_trim;                        % Subtract to find the sample length of all pulses

numBins = 10;
[N, edges, bin] = histcounts(enc_pulse_lengths,numBins);                    % Bin the pulse lengths (8 bins)
figure;histogram(enc_pulse_lengths,edges)
set(gca,'FontSize',12)
xlabel('Movement Duration (ms)','FontSize',14)
ylabel('Count','FontSize',14)
title('Movements binned by length','FontSize',14)
if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'Encoder_sequence_binned_durs.eps'))
end


preOnsetEnc = .3;     % Time before encoder onset to be incorporated (sec)
endTime = 2*1000;%max(edges);
timeWin = endTime+preOnsetEnc*1000;%max(edges)+preOnsetEnc*1000;                                      % Max(edges) determines the length of the window... this is = to the maximum histogram computed edge for sequence length

encoder_trials_binned = NaN(max(N),timeWin,length(N));
encChunkTime = -preOnsetEnc:0.001:1.999;%max(edges)/1000-.001;

disp('Plotting sorted encoder activity: Onsets')
numMvmt2plot = 4; %length(N)
co1 = othercolor('BuPu9',16);
co1 = co1(8:2:end,:);
numSubAxis = ceil(sqrt(numMvmt2plot));

mvmt_onset_sorted_fig = figure;
for oo = 1:numMvmt2plot %length(N)
    h3(oo) = subplot(numSubAxis,numSubAxis,oo); hold on
    if any(N(oo))
        temp = crossings2_trim(bin == oo);
        if any(temp > length(encoder_vel_norm_sub)-max(edges)-1)
            temp(temp > length(encoder_vel_norm_sub)-max(edges)-1) = [];
        end
        for jj = 1:length(temp)
            encoder_trials_binned(jj,1:length(encChunkTime),oo) = encoder_vel_norm_sub(temp(jj)-preOnsetEnc*1000: temp(jj)+ endTime-1);  % Dim 1: Trials,  Dim 2: time,  Dim 3: Encoder bout length
        end
        boundedline(encChunkTime,nanmean(encoder_trials_binned(:,:,oo)),nanstd(encoder_trials_binned(:,:,oo))/sqrt(N(oo)-1),'cmap',co1(oo,:),'transparency', 0.1,'alpha')
        set(gca,'FontSize',12)
    end
end
linkaxes(h3,'x');
suptitle('Encoder Velocity Sorted by Duration')

if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_Encoder_sequence_sorted_subplots_ON.eps'))
end

%% Repeat for aligning to the offset of movement
disp('Plotting sorted encoder activity: Offsets')

preOffsetEnc = 1;
encChunkTime_off = -preOffsetEnc:0.001:1-.001;
encoder_offset_trials_binned = NaN(max(N),length(encChunkTime_off),length(N));

mvmt_offset_sorted_fig = figure;
for oo2 = 1:numMvmt2plot
    h4(oo2) = subplot(numSubAxis,numSubAxis,oo2); hold on
    if any(N(oo2))
        temp2 = pulse_offsets(bin == oo2);
        if any(temp2 > length(encoder_vel_norm_sub)-max(edges)-1)
            temp2(temp2 > length(encoder_vel_norm_sub)-max(edges)-1) = [];
        end
        for jj = 1:length(temp2)
            encoder_offset_trials_binned(jj,1:length(encChunkTime_off),oo2) = encoder_vel_norm_sub(temp2(jj)-preOffsetEnc*1000: temp2(jj)+ 999);%max(edges)-1);  % Dim 1: Trials,  Dim 2: time,  Dim 3: Encoder bout length
        end
        boundedline(encChunkTime_off,nanmean(encoder_offset_trials_binned(:,:,oo2)),nanstd(encoder_offset_trials_binned(:,:,oo2))/sqrt(N(oo2)-1),'cmap',co1(oo2,:),'transparency', 0.1,'alpha')
        set(gca,'FontSize',12)
    end
end
linkaxes(h4,'x');
suptitle('Encoder Velocity Sorted by Duration')

if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_Encoder_sequence_sorted_subplots_OFF.eps'))
end


%% Do the same for EMG data
% Sort out different lengths of encoder activity
disp('Sorting EMG activity by duration')
timeEMG = linspace(0,length(EMG_smo(:,1))/1000,length(EMG_smo(:,1)));
fwd_window2 = 300;
EMG_thresh = [];
EMG_thresh = EMG_smo(:,1) > 5*bi_baseline_std;                      % Find where encoder is higher than threshold
ss = 1;
while ss <= length(EMG_smo(:,1)) - fwd_window2 
    if EMG_thresh(ss) == 1 && sum(EMG_thresh(ss:ss + fwd_window2)) > 1      % First combine encoder epochs that are within 'fwd_window2' of each other
        %if isempty(intersect(ss+1, crossings))                             % And don't overlap with detected crossings
            EMG_thresh(ss + 1) = 1;
        %end
    end
    ss = ss + 1;
end
if plotfigs
    figure; hold on % Check for accuracy
    plot(timeEMG,EMG_smo(:,1));
    scatter(timeEMG(crossings),EMG_smo(crossings,1))
    plot(timeEMG,EMG_thresh,'k','LineWidth',1.5)
    set(gca,'FontSize',12)
    xlabel('Time (sec)','FontSize',14)
    ylabel('EMG magnitude','FontSize',14)
    title('EMG sequence sorted','FontSize',14)
    if saveFigs
        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'EMG_sequence_chunking_examp.eps'))
    end

end
% For each encoder timestamp measure pulse duration; Bin & Histogram
pulse_onsets_EMG = [];
pulse_offsets_EMG=[];
pulse_onsets_EMG = find(diff(EMG_thresh) == 1);                                 % Find onsets
pulse_offsets_EMG = find(diff(EMG_thresh) == -1);                               % Find offsets

if pulse_offsets_EMG(1) < pulse_onsets_EMG(1)                                       % Trim pulse offsets that don't have an onset (onset before aquisition starts)
    pulse_offsets_EMG(1) = [];
end

if pulse_onsets_EMG(end) > pulse_offsets_EMG(end)                                   % Trim threshold crossing events w/o an offset (offset after aquisition ends)
    pulse_onsets_EMG(end) = [];
end

crossings_trim = pulse_onsets_EMG;
pulse_onsets_EMG_trim = [];
pulse_offsets_EMG_trim = [];
if length(crossings_trim) ~= length(pulse_offsets_EMG)
    for w = 1:length(pulse_onsets_EMG)
        if ~any(pulse_onsets_EMG(w) == crossings2)                              % If a crossing wasn't detected but there was a pulse
            pulse_onsets_EMG_trim = [pulse_onsets_EMG_trim; w];
            pulse_offsets_EMG_trim = [pulse_offsets_EMG_trim; w];
        end
    end
end
pulse_onsets_EMG(pulse_onsets_EMG_trim) = [];
pulse_offsets_EMG(pulse_offsets_EMG_trim) = [];

EMG_pulse_lengths = pulse_offsets_EMG - crossings_trim;                        % Subtract to find the sample length of all pulses

numBins = 10;
[N_EMG, edges_EMG, bin_EMG] = histcounts(EMG_pulse_lengths,numBins);                    % Bin the pulse lengths (8 bins)
figure;histogram(EMG_pulse_lengths,edges_EMG)
set(gca,'FontSize',12)
xlabel('Movement Duration (ms)','FontSize',14)
ylabel('Count','FontSize',14)
title('Movements binned by length','FontSize',14)
if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'EMG_sequence_binned_durs.eps'))
end


preOnsetEnc = .3;     % Time before  onset to be incorporated (sec)
timeWin_EMG = 2000+preOnsetEnc*1000;%max(edges_EMG)+preOnsetEnc*1000;                                      % Max(edges) determines the length of the window... this is = to the maximum histogram computed edge for sequence length

EMG_trials_binned = NaN(max(N_EMG),timeWin_EMG,length(N_EMG));
EMG_ChunkTime = -preOnsetEnc:0.001:1.999;%max(edges_EMG)/1000-.001;

disp('Plotting sorted EMG activity: Onsets')
numMvmt2plot = 4; %length(N)
co1 = othercolor('Greys8',16);
co1 = co1(8:2:end,:);
numSubAxis = ceil(sqrt(numMvmt2plot)); 
EMG_onset_sorted_fig = figure;
%get(mvmt_onset_sorted_fig)
for oo = 1:numMvmt2plot
    h3(oo) = subplot(numSubAxis,numSubAxis,oo); hold on
    %get(h3(oo));hold on;
    if any(N_EMG(oo)) && N_EMG(oo) > 1
        temp = crossings_trim(bin_EMG == oo);
        if any(temp > length(EMG_smo(:,1))-max(edges_EMG)-1)
            temp(temp > length(EMG_smo(:,1))-max(edges_EMG)-1) = [];
        end
        for jj = 1:length(temp)
            EMG_trials_binned(jj,1:length(EMG_ChunkTime),oo) = EMG_smo(temp(jj)-preOnsetEnc*1000: temp(jj)+ endTime-1,1);  % Dim 1: Trials,  Dim 2: time,  Dim 3: EMG bout length
        end
        boundedline(EMG_ChunkTime,nanmean(EMG_trials_binned(:,:,oo)),nanstd(EMG_trials_binned(:,:,oo))/sqrt(N_EMG(oo)-1),'cmap',co1(oo,:),'transparency', 0.1,'alpha')
        set(gca,'FontSize',12)
    end
    hold off;
end
linkaxes(h3,'x');
suptitle('EMG Activity Sorted by Duration')

if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_EMG_sequence_sorted_subplots_ON.eps'))
end
%%
%%%%                                                                 %%%%
%                                                                       %
%   Populate 4d array w/ Ca data relative to encoder movement ONSET     %
%                                                                       %
%%%%                                                                 %%%% 

% Dimension 1: Trial #
% Dimension 2: Time
% Dimension 3: Cell #
% Dimension 4: Sequence length

disp('Aligning Ca data to sorted encoder data')
encoder_Ca_seq_sorted = NaN(max(N),max(edges)+preOnsetEnc*1000,size(pruned_Cell_ts_ca,2),length(N));
encoder_Sp_seq_sorted = NaN(max(N),max(edges)+preOnsetEnc*1000,size(pruned_Cell_ts_sp,2),length(N));

for aa = 1:length(N)                                                        % length of N corresponds to the number of different movement epochs (defined by their length)
    current_trials = crossings2_trim(bin == aa);                            % Indices of trial onset for current sequence length (the indices of encoder crossings for each trial of a given bin of movement length)
    
    if any(current_trials > length(encoder_vel_norm_sub)-max(edges)-1)      % Get rid of trials too close to end of acquisition
       current_trials(current_trials > length(encoder_vel_norm_sub)-max(edges)-1) = [];
    end
    
    for bb = 1:length(current_trials)                                       % Number of trials available for that sequence length
        for cc = 1:size(pruned_Cell_ts_ca,2)
            encoder_Ca_seq_sorted(bb,1:max(edges)+preOnsetEnc*1000,cc,aa) = pruned_Cell_ts_ca(current_trials(bb)-preOnsetEnc*1000 : current_trials(bb)+max(edges)-1,cc); % Ca activity for each ROI, for each movement, for each movement length bin
            encoder_Sp_seq_sorted(bb,1:max(edges)+preOnsetEnc*1000,cc,aa) = pruned_Cell_ts_sp(current_trials(bb)-preOnsetEnc*1000 : current_trials(bb)+max(edges)-1,cc); % Sp activity for each ROI, for each movement, for each movement length bin        
        end
    end
    
end

% Conditions to specify for plotting
numCellsToPlot = size(pruned_Cell_ts_ca,2);
numMvmtsToPlot = 3;                                                         % The movments to plot
encChunkTime_on = -preOnsetEnc:0.001:edges(numMvmtsToPlot+1)/1000-.001;     % Remake time so extra crap isn't added
numSubAxis = ceil(sqrt(numCellsToPlot));

% Plot each cell's mean Ca activity relative to sequence length, ONSET
%co1 = othercolor('RdPu9',16);
% co1 = othercolor('YlOrBr9',18);
% co1 = co1(8:3:end,:);
co1 = othercolor('BuPu9',18);
co1 = co1(8:3:end,:);
%co1 = co1(6:5:end,:);
if plotfigs
    spikes = 0;
    norm = 0;
    figure;
    for aa = 1:numCellsToPlot
        h4(aa) = subplot(numSubAxis,numSubAxis,aa);
        for bb = 1:numMvmtsToPlot                                           % Only subplot the first bb sequence lengths
            
            %aa = i_order(aa,1);
            
            current_sequence_ca_on = encoder_Ca_seq_sorted(:,1:length(encChunkTime_on),:,bb);
            current_sequence_sp_on = encoder_Sp_seq_sorted(:,1:length(encChunkTime_on),:,bb);
            
            if ~spikes && ~norm     % Ca activity
                boundedline(encChunkTime_on, nanmean(current_sequence_ca_on(:,:,aa)),nanstd(current_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'Ca Activity Relative to Movement Duration, ONSET (colors)';
            elseif spikes           % Spikes
                boundedline(encChunkTime_on, nanmean(current_sequence_sp_on(:,:,aa)),nanstd(current_sequence_sp_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'Sp Activity Relative to Movement Duration, ONSET (colors)';
            elseif norm             % Feature scaled Ca activity
                boundedline(encChunkTime_on, (nanmean(current_sequence_ca_on(:,:,aa))-min(nanmean(current_sequence_ca_on(:,:,aa))))/(max(nanmean(current_sequence_ca_on(:,:,aa)))-min(nanmean(current_sequence_ca_on(:,:,aa)))),...
                    nanstd(current_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
%                 boundedline(encChunkTime_on, (nanmean(current_sequence_ca_on(:,:,aa))-min(nanmean(current_sequence_ca_on(:,:,aa))))/(max(nanmean(current_sequence_ca_on(:,:,aa)))-min(nanmean(current_sequence_ca_on(:,:,aa)))),...
%                     nanstd(current_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
%                 normCa_curr_cell = [];
%                 normCa_curr_cell = (current_sequence_ca_on(:,:,aa)-min(current_sequence_ca_on(:,:,aa)))/((max(current_sequence_ca_on(:,:,aa))-min(current_sequence_ca_on(:,:,aa));
%                 
                title4 = 'Scaled ROI Ca Activity Relative to Movement Duration, ONSET (colors),SEM';
            end
        end
    end
    suptitle(title4)
    
    if saveFigs
        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_ROI_Ca_Act_Relative_Mvmt_Duration_ON.eps'))
    end
    linkaxes(h4,'x')
end

%%
%%%%                                                        %%%%
%                                                              %
%   Populate 4d array w/ Ca data relative to EMG ONSET    %
%                                                              %
%%%%                                                        %%%% 

% Dimension 1: Trial #
% Dimension 2: Time
% Dimension 3: Cell #
% Dimension 4: Sequence length

disp('Aligning Ca data to sorted encoder data')
EMG_Ca_seq_sorted = NaN(max(N_EMG),timeWin_EMG,size(pruned_Cell_ts_ca,2),length(N_EMG));
EMG_Sp_seq_sorted = NaN(max(N_EMG),timeWin_EMG,size(pruned_Cell_ts_sp,2),length(N_EMG));
current_EMG_trials=[];
for aa = 1:length(N_EMG)                                                        % length of N corresponds to the number of different movement epochs (defined by their length)
    current_EMG_trials = crossings_trim(bin_EMG == aa);                            % Indices of trial onset for current sequence length (the indices of encoder crossings for each trial of a given bin of movement length)
    
    if any(current_EMG_trials > length(EMG_smo(:,1))-max(edges_EMG)-1)      % Get rid of trials too close to end of acquisition
       current_EMG_trials(current_EMG_trials > length(EMG_smo(:,1))-max(edges_EMG)-1) = [];
    end
    
    for bb = 1:length(current_EMG_trials)                                       % Number of trials available for that sequence length
        for cc = 1:size(pruned_Cell_ts_ca,2)
            EMG_Ca_seq_sorted(bb,1:timeWin_EMG,cc,aa) = pruned_Cell_ts_ca(current_EMG_trials(bb)-preOnsetEnc*1000 : current_EMG_trials(bb)+max(edges_EMG)-1,cc); % Ca activity for each ROI, for each movement, for each movement length bin
            EMG_Sp_seq_sorted(bb,1:timeWin_EMG,cc,aa) = pruned_Cell_ts_sp(current_EMG_trials(bb)-preOnsetEnc*1000 : current_EMG_trials(bb)+max(edges_EMG)-1,cc); % Sp activity for each ROI, for each movement, for each movement length bin        
        end
    end
    
end

% Conditions to specify for plotting
numCellsToPlot = size(pruned_Cell_ts_ca,2);
numMvmtsToPlot = 3;                                                         % The movments to plot
EMGChunkTime_on = -preOnsetEnc:0.001:edges_EMG(numMvmtsToPlot+1)/1000-.001;     % Remake time so extra crap isn't added
numSubAxis = ceil(sqrt(numCellsToPlot));

% Plot each cell's mean Ca activity relative to sequence length, ONSET
%co1 = othercolor('RdPu9',16);
% co1 = othercolor('YlOrBr9',18);
% co1 = co1(8:3:end,:);
co1 = othercolor('BuPu9',18);
co1 = co1(8:3:end,:);
%co1 = co1(6:5:end,:);
if plotfigs
    spikes = 0;
    norm = 0;
    figure;
    for aa = 1:numCellsToPlot
        h4(aa) = subplot(numSubAxis,numSubAxis,aa);
        for bb = 1:numMvmtsToPlot                                           % Only subplot the first bb sequence lengths
            
            %aa = i_order(aa,1);
            
            current_EMG_sequence_ca_on = EMG_Ca_seq_sorted(:,1:length(EMGChunkTime_on),:,bb);
            current_EMG_sequence_sp_on = EMG_Sp_seq_sorted(:,1:length(EMGChunkTime_on),:,bb);
            
            if ~spikes && ~norm     % Ca activity
                boundedline(EMGChunkTime_on, nanmean(current_EMG_sequence_ca_on(:,:,aa)),nanstd(current_EMG_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'Ca Activity Relative to Movement Duration, ONSET (colors)';
            elseif spikes           % Spikes
                boundedline(EMGChunkTime_on, nanmean(current_EMG_sequence_sp_on(:,:,aa)),nanstd(current_EMG_sequence_sp_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'Sp Activity Relative to Movement Duration, ONSET (colors)';
            elseif norm             % Feature scaled Ca activity
                boundedline(EMGChunkTime_on, (nanmean(current_EMG_sequence_ca_on(:,:,aa))-min(nanmean(current_EMG_sequence_ca_on(:,:,aa))))/(max(nanmean(current_EMG_sequence_ca_on(:,:,aa)))-min(nanmean(current_EMG_sequence_ca_on(:,:,aa)))),...
                    nanstd(current_EMG_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
%                 boundedline(encChunkTime_on, (nanmean(current_sequence_ca_on(:,:,aa))-min(nanmean(current_sequence_ca_on(:,:,aa))))/(max(nanmean(current_sequence_ca_on(:,:,aa)))-min(nanmean(current_sequence_ca_on(:,:,aa)))),...
%                     nanstd(current_sequence_ca_on(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
%                 normCa_curr_cell = [];
%                 normCa_curr_cell = (current_sequence_ca_on(:,:,aa)-min(current_sequence_ca_on(:,:,aa)))/((max(current_sequence_ca_on(:,:,aa))-min(current_sequence_ca_on(:,:,aa));
%                 
                title4 = 'Scaled ROI Ca Activity Relative to EMG Duration, ONSET (colors),SEM';
            end
        end
    end
    suptitle(title4)
    
    if saveFigs
        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_ROI_Ca_Act_Relative_EMG_Duration_ON.eps'))
    end
    linkaxes(h4,'x')
end
%%
%%%%                                                        %%%%
%                                                              %
%   Populate 4d array w/ Ca data relative to movement OFFSET   %
%                                                              %
%%%%                                                        %%%% 


disp('Aligning Ca data to sorted encoder data: offset')
encoder_Ca_seq_sorted_offset = NaN(max(N),length(encChunkTime_off),size(pruned_Cell_ts_ca,2),length(N));
encoder_Sp_seq_sorted_offset = NaN(max(N),length(encChunkTime_off),size(pruned_Cell_ts_sp,2),length(N));

for aa = 1:length(N)                                                        % length of N corresponds to the number of different movement epochs (defined by their length)
    current_trials = pulse_offsets(bin == aa);                              % Indices of trial onset for current sequence length (the indices of encoder crossings for each trial of a given bin of movement length)
    
    if any(current_trials > length(encoder_vel_norm_sub)-max(edges)-1)
       current_trials(current_trials > length(encoder_vel_norm_sub)-max(edges)-1) = [];
    end
    
    for bb = 1:length(current_trials)                                       % Number of trials available for that sequence length
        for cc = 1:size(pruned_Cell_ts_ca,2)
            encoder_Ca_seq_sorted_offset(bb,1:length(encChunkTime_off),cc,aa) = pruned_Cell_ts_ca(current_trials(bb)-preOffsetEnc*1000 : current_trials(bb)+999,cc); % Ca activity for each ROI, for each movement, for each movement length bin
            encoder_Sp_seq_sorted_offset(bb,1:length(encChunkTime_off),cc,aa) = pruned_Cell_ts_sp(current_trials(bb)-preOffsetEnc*1000 : current_trials(bb)+999,cc); % Sp activity for each ROI, for each movement, for each movement length bin        
        end
    end
    
end

% Conditions to specify for plotting
numCellsToPlot = size(pruned_Cell_ts_ca,2);
numMvmtsToPlot = 3;                                                         % The movments to plot
encChunkTime_off = -preOffsetEnc:0.001:.999;                                % Remake time so extra crap isn't added
numSubAxis = ceil(sqrt(numCellsToPlot));

% Plot each cell's mean Ca activity relative to sequence length, OFFSET
if plotfigs
    figure;
    spikes = 0;
    norm = 1;
    for aa = 1:numCellsToPlot
        h4(aa) = subplot(numSubAxis,numSubAxis,aa);
        for bb = 1:numMvmtsToPlot                                           % Only subplot the first bb sequence lengths
            
            current_sequence_ca_off = encoder_Ca_seq_sorted_offset(:,1:length(encChunkTime_off),:,bb);
            current_sequence_sp_off = encoder_Sp_seq_sorted_offset(:,1:length(encChunkTime_off),:,bb);
            
            if ~spikes && ~norm         % Ca activity
                boundedline(encChunkTime_off, nanmean(current_sequence_ca_off(:,:,aa)),nanstd(current_sequence_ca_off(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'OFFSET Aligned Ca Activity';
            elseif spikes               % Spikes
                boundedline(encChunkTime_off, nanmean(current_sequence_sp_off(:,:,aa)),nanstd(current_sequence_sp_off(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'OFFSET Aligned Sp Activity';
            elseif norm                 % Feature scaled Ca activity
                boundedline(encChunkTime_off, (nanmean(current_sequence_ca_off(:,:,aa))-min(nanmean(current_sequence_ca_off(:,:,aa))))/(max(nanmean(current_sequence_ca_off(:,:,aa)))-min(nanmean(current_sequence_ca_off(:,:,aa)))),...
                    nanstd(current_sequence_ca_off(:,:,aa))/sqrt(N(bb)-1),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
                title4 = 'Scaled OFFSET Aligned Ca Activity, SEM';
            end
            
        end
    end
    suptitle(title4)
    
    if saveFigs
        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_ROI_Ca_Act_Relative_Mvmt_Duration_OFF.eps'))
    end
    linkaxes(h4,'x')
end

%% Sort and heatmap normalized Ca activity

% Normalize matrix of activity so each cell is 0-1






%%% Looks like longer movements peak in amplitude earlier, and the Ca
%%% activity relative to thos
encoder_Ca_seq_sorted_onset_squeezed = squeeze(nanmean(encoder_Ca_seq_sorted(:,1:1500,:,1:4),1));
    encoder_Ca_seq_sorted_onset_min = repmat(nanmin(encoder_Ca_seq_sorted_onset_squeezed),size(encoder_Ca_seq_sorted_onset_squeezed,1),1,1);
    encoder_Ca_seq_sorted_onset_max = repmat(nanmax(encoder_Ca_seq_sorted_onset_squeezed),size(encoder_Ca_seq_sorted_onset_squeezed,1),1,1);
    encoder_Ca_seq_sorted_onset_norm = (encoder_Ca_seq_sorted_onset_squeezed - encoder_Ca_seq_sorted_onset_min)./(encoder_Ca_seq_sorted_onset_max - encoder_Ca_seq_sorted_onset_min);
% Average across movement conditions (normalized), and heatmap

norm_Ca_onset_reorder = nan(size(encoder_Ca_seq_sorted_onset_squeezed));
for i = 1:size(encoder_Ca_seq_sorted_onset_norm,3)
    [~,i_max] = max(encoder_Ca_seq_sorted_onset_norm(:,:,i));
    [~, i_order(:,i)] = sort(i_max);
    norm_Ca_onset_reorder(:,:,i) = encoder_Ca_seq_sorted_onset_norm(:,i_order(:,i),i);
end
figure; imagesc(norm_Ca_onset_reorder(:,:,1)')
figure; imagesc(norm_Ca_onset_reorder(:,:,2)')
figure; imagesc(norm_Ca_onset_reorder(:,:,3)')
figure; imagesc(norm_Ca_onset_reorder(:,:,4)')

        cd(figDir)
        print('-painters', '-depsc', strcat(CNMF_filename(9:23),'_Heatmap_Ca_Activity_Movement_ON.eps'))

mean_norm_Ca_onset = mean(encoder_Ca_seq_sorted_onset_norm,3)';
[~,i_max] = max(mean_norm_Ca_onset,[],2);
[~, i_order_m] = sort(i_max);
mean_norm_Ca_onset_reorder = mean_norm_Ca_onset(i_order_m,:);
figure;
imagesc(mean_norm_Ca_onset_reorder)


%% Save

analyzed.encoder_trials_binned = encoder_trials_binned;
analyzed.encoder_offset_trials_binned = encoder_offset_trials_binned;
analyzed.encChunkTime_on = encChunkTime_on;
analyzed.encChunkTime_off = encChunkTime_off;
analyzed.numCellsToPlot = numCellsToPlot;
analyzed.numMvmtsToPlot = numMvmtsToPlot;                                                                      
analyzed.numSubAxis = numSubAxis;
analyzed.encoder_Ca_seq_sorted_offset = encoder_Ca_seq_sorted_offset;
analyzed.encoder_Sp_seq_sorted_offset = encoder_Sp_seq_sorted_offset;
analyzed.encoder_Ca_seq_sorted = encoder_Ca_seq_sorted;
analyzed.encoder_Sp_seq_sorted = encoder_Sp_seq_sorted;

analyzed.encoder_Ca_seq_sorted_onset_norm = encoder_Ca_seq_sorted_onset_norm;

cd(analyzDir)
saveName = strcat(CNMF_filename(1:end-4), '_analyzed');
save(saveName,'-struct','analyzed')

%%
%%%%     %%%%
%           %
%   Extra   %
%           %
%%%%     %%%% 
%% Make movie of encoder, etc data relative to movie playback

load(fullfile(mov_dir,mov_filename));
Y = mov_nonrigid_corr;
[d1,d2,T] = size(Y);     
d = d1*d2;                                         
Yr = reshape(Y,d,T);

timeV = linspace(0,length(encoder)/1000,length(encoder_vel_norm_sub));
encoder_vel_norm_sub_pad = [zeros(1000,1); encoder_vel_norm_sub; zeros(6000,1)];
Cell_ts_ca_pad = [zeros(1000,size(Cell_ts_ca,2)); Cell_ts_ca; zeros(6000,size(Cell_ts_ca,2))];
fr = (length(encoder_vel_norm_sub)/1000);
sr = length(f)/fr;
make_patch_video_edit(A2,C2,b2,f2,Yr,Coor,options)%,encoder_vel_norm_sub_pad,Cell_ts_ca_pad,sr)
%% Run cross correlation between Ca activity and movement
% A = 4; % Cell
% B = 1; % Movement length
% 
% activ = nanmean(encoder_Ca_seq_sorted(:,:,A,B));
% mvmt = nanmean(encoder_trials_binned(:,:,B));
% figure; hold on;
% plot(activ,'b')
% plot(mvmt,'r')

% Try this: correlate Ca activity to 'non-corresponding' movement epoch
% lengths
 
% Need to preallocate arrays with nans.. without doing this the means will
% be screwed up
maxLag = .5;
xcorr_trials=[]; %= nan(size(encoder_Ca_seq_sorted,1),size(encoder_Ca_seq_sorted,2)*2-1,size(encoder_Ca_seq_sorted,3),size(encoder_Ca_seq_sorted,4));
xcorr_trials_lag=[]; %= nan(size(encoder_Ca_seq_sorted,1),size(encoder_Ca_seq_sorted,2)*2-1,size(encoder_Ca_seq_sorted,3),size(encoder_Ca_seq_sorted,4));
xcorr_peaks = nan(size(encoder_Ca_seq_sorted,1),size(encoder_Ca_seq_sorted,3),size(encoder_Ca_seq_sorted,4));
xcorr_peaks_idx = nan(size(encoder_Ca_seq_sorted,1),size(encoder_Ca_seq_sorted,3),size(encoder_Ca_seq_sorted,4));

for aa = 1:length(N) % length of N corresponds to the number of different movement epochs (defined by their length)
    current_trials = crossings2_trim(bin == aa);     % Indices of trial onset for current sequence length (the indices of encoder crossings for each trial of a given bin of movement length)
    for bb = 1:length(current_trials)    % Number of trials available for that sequence length
        for cc = 1:size(pruned_Cell_ts_ca,2)
            [xcorr_trials(bb,:,cc,aa), lag] = xcorr(encoder_Ca_seq_sorted(bb,:,cc,aa),nanmean(encoder_trials_binned(:,:,aa)),maxLag*1000); % Ca activity precedes the encoder when the peak is < 0
            [xcorr_peaks(bb,cc,aa), xcorr_peaks_idx(bb,cc,aa)] = max(xcorr_trials(bb,:,cc,aa),[],2);
        end
    end
end
% Make time for cross correlation lag
time_xcorr = linspace(-length(xcorr_trials)/1000/2,length(xcorr_trials)/1000/2,length(xcorr_trials));


% Try to get mean lag info
xcorr_peaks_lag_mean = squeeze(nanmean(xcorr_peaks_idx,1));
xcorr_peaks_lag_std = squeeze(nanstd(xcorr_peaks_idx,1));

figure; hold on;
for x = 1:3
    errorbar(xcorr_peaks_lag_mean(:,x),xcorr_peaks_lag_std(:,x)/sqrt(size(xcorr_peaks_idx,1)));
end
    plot([0 size(xcorr_peaks_idx,2)], [find(lag==0) find(lag==0)]);
        
numCellsToPlot = size(pruned_Cell_ts_ca,2);
numSubAxis = ceil(sqrt(numCellsToPlot));
figure;
for aa = 1:numCellsToPlot
    h5(aa) = subtightplot(numSubAxis,numSubAxis,aa,[.02,.02]);
    for bb = 1:3 % Only subplot the first 4 sequence lengths
        current_sequence_xcorr = xcorr_trials(:,:,:,bb);
        %plot(lag,nanmean(current_sequence_xcorr(:,:,aa))/max(nanmean(current_sequence_xcorr(:,:,aa)))); hold on;
        plot(lag,nanmean(current_sequence_xcorr(:,:,aa))); hold on;
        %boundedline(lag, nanmean(current_sequence_xcorr(:,:,aa))/max(nanmean(current_sequence_xcorr(:,:,aa))),nanstd(current_sequence_xcorr(:,:,aa))/sqrt(size(current_sequence_xcorr,1))/max(nanmean(current_sequence_xcorr(:,:,aa))),'cmap',co1(bb,:),'transparency', 0.1,'alpha'); hold on
    end
    plot([0 0],[0 1.5])
end
suptitle('Normalized ROI v. Movement Cross Correlation Relative to Movement Duration (colors)')
if saveFigs
    cd(figDir)
    print -painters -depsc Encoder_sequence_sorted_Ca.eps 
end
linkaxes(h5,'x')
linkaxes(h5,'y')
%% Identify rabies-labeled ROIs
    if any(rabies_img)
        options.rabies_thresh = 1;
        rabies_red = imread(fullfile(rabies_red_dir,rabies_img));
    end
figure;
[Coor,json_file] = plot_contours(A2_pruned,rabies_red,options,1);

% Generate mask
rabies_red_mask = rabies_red > 100; % Generate mask from rabies labeling.  Set threshold for keeping. 150 initially
figure; imagesc(rabies_red_mask)
xx = double(medfilt2(rabies_red_mask));
xx(xx == 0) = nan;
figure;imagesc(xx)

rabies_img_med = double(medfilt2(rabies_red));
rabies_red_mask2 = rabies_img_med >150;
figure;imagesc(rabies_red_mask2);


figure;
[Coor,json_file] = plot_contours(A2_pruned,xx,options,1);

% Mask pixels that correspond to each ROI
ROImask = cell(1,length(Coor));
ROImaskSum = nan(256,256);
for o = 1:length(Coor)
    ROImask{o} = double(poly2mask(Coor{o,:}(1,:),Coor{o,:}(2,:),256,256));
    ROImask{o}(ROImask{o} == 0) = nan;
    ROImaskSum = ROImaskSum + ROImask{o}; % To visualize all masks collapsed in 1 image
end
figure; imagesc(ROImask{1})

colabeled_pix = cell(1,length(Coor));
colabel_thresh = 60; % Set percent threshold for defining substantial colabeling of ROI & rabies
ROI_dbl_labeled = [];
% Run through each mask to see if there are rabies pixels within it
for oo = 1:length(Coor)
     colabeled_pix{oo} = find((ROImask{1,oo} - rabies_red_mask) == 0);
    if ((length(colabeled_pix{1,oo}) - 1) / nansum(nansum(ROImask{oo}))) >= colabel_thresh 
        ROI_dbl_labeled = [ROI_dbl_labeled oo];
    end
end





% %%%%%% Some extra crap
% %% Select ROIs by their average change in activity
% disp('Sorting ROIs by their activity')
% mean_bi_crossings = mean(bi_crossings_Ca);
% keepers_inc = [];
% for ii = 1:size(mean_bi_crossings,3)
%     if mean(mean_bi_crossings(1,pre_win_init:end,ii)) > 5*mean(mean_bi_crossings(1,1:pre_win_init,ii))
%         keepers_inc = [keepers_inc; ii];
%     end
% end
%     bi_crossings_inc = [];
%     bi_Ca_incs = bi_crossings_Ca(:,:,keepers_inc);
%     
% numSubAxis = ceil(sqrt(length(keepers_inc)));
% figure; 
% for cc = 1:length(keepers_inc)
%     h6(cc) = subplot(numSubAxis,numSubAxis,cc); hold on
%     boundedline(time_init,mean(bi_Ca_incs(:,:,cc))/max(mean(bi_Ca_incs(:,:,cc))),std(bi_Ca_incs(:,:,cc)/max(mean(bi_Ca_incs(:,:,cc))))/sqrt(length(keepers_inc)),'-b','transparency', 0.03,'alpha')
%     ylabel('Ca activity [norm to max]','fontsize',9)
%     xlabel('Time relative to onset (msec)','fontsize',9)
% end
% linkaxes(h6,'x')
% suptitle('Selected ROIs: Ca data aligned to bi EMG onset')
% 
% % Select ROIs by their average change in Z score before/after EMG onset
% keepers_inc = [];
% for ii = 1:size(mu_bi_Ca,3)
%     if mean(mu_bi_Ca(1,pre_win_init:end,ii)) > 5*mean(mu_bi_Ca(1,1:pre_win_init,ii))
%         keepers_inc = [keepers_inc; ii];
%     end
% end
%     mu_bi_Ca_incs = [];
%     sig_bi_Ca_incs = [];
%     mu_bi_Ca_incs = mu_bi_Ca(:,:,keepers_inc);
%     sig_bi_Ca_incs = sig_bi_Ca(:,:,keepers_inc);
% keepers_dec = [];
% for ii = 1:size(mu_bi_Ca,3)
%     if mean(mu_bi_Ca(1,1:pre_win_init,ii)) > 1.5*mean(mu_bi_Ca(1,pre_win_init:end,ii))
%         keepers_dec = [keepers_dec; ii];
%     end
% end
%     mu_bi_Ca_dec = [];
%     sig_bi_Ca_dec = [];
%     mu_bi_Ca_dec = mu_bi_Ca(:,:,keepers_dec);
%     sig_bi_Ca_dec = sig_bi_Ca(:,:,keepers_dec);    
%       
% numSubAxis = ceil(sqrt(length(keepers_inc)));
% figure; 
% for cc = 1:length(keepers_inc)
%     h(cc) = subplot(numSubAxis,numSubAxis,cc); hold on
%     boundedline(time_init,mu_bi_Ca_incs(:,:,cc),sig_bi_Ca_incs(:,:,cc)/sqrt(length(keepers_inc)),'-b','transparency', 0.03,'alpha')
%     ylabel('Ca activity [Z Score]','fontsize',9)
%     xlabel('Time relative to onset (msec)','fontsize',9)
% end
% linkaxes(h,'x')
% suptitle('Selected ROIs: Ca data aligned to bi EMG onset')    
% 
% %% Align Data by Reward
% % Align data by pulse train onset, organized by # of pulses
% uniqPulseTrains = unique(rew_info);
% uniqPulseTrains(uniqPulseTrains == 0) = [];
% 
% pre_onset = 1; % Windows for viewing data relative to reward onset
% post_onset = 1;
% 
% % Remove pulse trains close to session on or off from analysis
% %%%% Stupid way of doing this.... rew_on < pre_onset = [];
% pulseTrains_valid = rew_info;
% if any(find(rew_on) < pre_onset*samp_rate_dwnsmp)
%     pulseTrains_valid(find(rew_on) < pre_onset*samp_rate_dwnsmp) = [];      
% end
% if any(find(rew_on) > length(rew_on) - post_onset*samp_rate_dwnsmp)
%     pulseTrains_valid(find(rew_on) > length(rew_on) - post_onset*samp_rate_dwnsmp) = [];
% end
% 
% 
% % Get onset indices of each reward pulse train length (1-5)
% rew_train_indices = cell(1,5);
% Ca_rew_train = cell(1,5);
%     clear Ca_rew_train_means
%     Ca_rew_train_means = cell(1,5);
%     clear Ca_rew_train_std
%     Ca_rew_train_std = cell(1,5);
% for pp = 1:length(uniqPulseTrains) 
%     rew_train_indices{uniqPulseTrains(pp)} = find(rew_info == uniqPulseTrains(pp));
%     
%     % Populate array with Data corresponding to each reward event, timed around pre/post onset specification
%     for jj = 1:length(rew_train_indices{uniqPulseTrains(pp)}) % The number of each rew train type
%             Ca_rew_train{1,uniqPulseTrains(pp),jj} = Cell_ts_ca(rew_train_indices{1,uniqPulseTrains(pp)}(jj) - pre_onset*samp_rate_dwnsmp : ...
%                 rew_train_indices{1,uniqPulseTrains(pp)}(jj) + post_onset*samp_rate_dwnsmp,:);
%     end
%      Ca_rew_train_means{1,uniqPulseTrains(pp)} = mean(cell2mat(Ca_rew_train(1,uniqPulseTrains(pp),:)),3);
%      Ca_rew_train_std{1,uniqPulseTrains(pp)} = std(cell2mat(Ca_rew_train(1,uniqPulseTrains(pp),:)),0,3);
% 
% end

EMG_smo_bi_norm = (EMG_smo(:,1) - min(EMG_smo(:,1)) )./ (max(EMG_smo(:,1)) - min(EMG_smo(:,1)));
%%Plot Encoder and EMG sorted segments on top of each other
disp('Plotting sorted EMG And Ecoder activity: Onsets')
numMvmt2plot = 4; %length(N)
co1 = othercolor('BuPu9',16);
co1 = co1(8:2:end,:);
co2 = othercolor('Greys8',16);
co2 = co2(8:2:end,:);

numSubAxis = ceil(sqrt(numMvmt2plot));
EMG_Encoder_onset_sorted_fig = figure;
%get(mvmt_onset_sorted_fig)
for oo = 1:numMvmt2plot
    h3(oo) = subplot(numSubAxis,numSubAxis,oo); hold on
    %get(h3(oo));hold on;
    if any(N_EMG(oo)) && N_EMG(oo) > 1
        temp = crossings_trim(bin_EMG == oo);
        if any(temp > length(EMG_smo(:,1))-max(edges_EMG)-1)
            temp(temp > length(EMG_smo(:,1))-max(edges_EMG)-1) = [];
        end
        for jj = 1:length(temp)
            encoder_trials_binned(jj,1:length(encChunkTime),oo) = encoder_vel_norm_sub(temp(jj)-preOnsetEnc*1000: temp(jj)+ endTime-1);  % Dim 1: Trials,  Dim 2: time,  Dim 3: Encoder bout length
            EMG_trials_binned(jj,1:length(EMG_ChunkTime),oo) = EMG_smo_bi_norm(temp(jj)-preOnsetEnc*1000: temp(jj)+ endTime-1);  % Dim 1: Trials,  Dim 2: time,  Dim 3: EMG bout length
        end
        boundedline(encChunkTime,nanmean(encoder_trials_binned(:,:,oo)),nanstd(encoder_trials_binned(:,:,oo))/sqrt(N(oo)-1),'cmap',co1(oo,:),'transparency', 0.1,'alpha')
        boundedline(EMG_ChunkTime,nanmean(EMG_trials_binned(:,:,oo)),nanstd(EMG_trials_binned(:,:,oo))/sqrt(N_EMG(oo)-1),'cmap',co2(oo,:),'transparency', 0.1,'alpha')
        set(gca,'FontSize',12)
        set(gca,'ylim',[0 0.3])
    end
    hold off;
end
linkaxes(h3,'x');
suptitle('Encoder(color) EMG (grey) Activity Sorted by Duration')
%%

if saveFigs
    cd(figDir)
    print('-painters', '-depsc', strcat(CNMF_filename(9:23),'Enc_EMG_sorted_onset.eps'))
end

end




























%%%%%%%%%%%%%%%%%%%%%%%


% % Feature scale mean neuron activity
% mean_crossings_Ca = mean(bi_crossings_Ca);
% mean_crossings_Ca_norm = (mean_crossings_Ca - min(mean_crossings_Ca,[],2)) ./ repmat((max(mean_crossings_Ca,[],2) - min(mean_crossings_Ca,[],2)),1,length(mean_crossings_Ca));     
% SE_crossings_Ca_norm = std(bi_crossings_Ca)/sqrt(num_init_events);

% % Normalize mean neuron activity
% mean_enc_crossings_Ca = mean(enc_crossings_Ca);
% mean_crossings_Ca2_norm = (mean_enc_crossings_Ca - min(mean_enc_crossings_Ca,[],2)) ./ repmat((max(mean_enc_crossings_Ca,[],2) - min(mean_enc_crossings_Ca,[],2)),1,length(mean_crossings_Ca));     
% SE_crossings_Ca2_norm = std(enc_crossings_Ca)/sqrt(num_init_events);

