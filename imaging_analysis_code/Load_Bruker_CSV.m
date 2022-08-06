%%%% Load CSV files collected on Bruker machine
data_path = 'G:\CSMNs\D12';
csv_filename = 'TSeries-09282017-D12-001_Cycle00001_VoltageRecording_001.csv';
CNMF_filename = 'D12_09282017_001_CNMF_results.mat';

%Set Params
samp_rate = 10000; %original sampling rate. Probably 10000
Y_galvo_col = 2; %column in csv of Y_galvo signal
encoder_col = 3; %column in csv of encoder signal
solenoid_col = 4; %column in csv of solenoid signal
EMG_cols = 5:8; %columns of EMG signals in csv
lp = 40; %low pass filter freq (Hz) for EMG filtering
smoFiltSTD = 10; %STD of Gaussian for EMG filtering
sub_rate = 10; %rate at which to downsample data - to go from 10K to 1k, use 10
Tri_ch = 2; %EMG header position where the triceps signal comes from
Bi_ch = 1; %EMG header position where the biceps signal comes from
event_dur = 0.5; % duration of alternation event epochs, in sec
modal_event_set_size = 30;

%Read CSV
cd(data_path)
csv_readtime = tic;
Mat = csvread(csv_filename,1,0) ;
% M = csvread(filename,R1,C1,[R1 C1 R2 C2]) 
toc(csv_readtime)
% 51 seconds for 20 minute movie

%Downsample and separate signals
Y_galvo = downsample(Mat(:,Y_galvo_col),sub_rate);
encoder = downsample(Mat(:,encoder_col),sub_rate);
solenoid = downsample(Mat(:,solenoid_col),sub_rate);
EMG = downsample(Mat(:,EMG_cols),sub_rate);

%Filter the EMG signals
EMG_smo = 0*EMG;
for i = 1:size(EMG,2)
    EMG_ch_time = tic;
	EMG_smo(:,i) = filterEMG(EMG(:,i),lp,smoFiltSTD);
    toc(EMG_ch_time)
end

%Subtract the modal value off EMG signals to make baseline 0 - use last 5 minutes, when there is often not much activity
for i = 1:size(EMG,2)
   num_datapoints = 300*samp_rate/sub_rate;
   [x,n] = hist(EMG_smo(end-num_datapoints:end,i),500);
   [m,ind] = max(x);
   signalMode = n(ind);
   EMG_smo(:,i) = EMG_smo(:,i) - signalMode;
end

%Normalize the EMG values to their 99.5th% value
for i = 1:size(EMG,2)
   m = quantile(EMG_smo(:,i),0.995);
   EMG_smo(:,i) = EMG_smo(:,i)/m;
   figure; plot(EMG_smo(:,i))
end

%Find alternation events
pre_win = floor(event_dur*1000/2);
[~, locs_tri] = findpeaks(EMG_smo(:,Tri_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
locs_tri(locs_tri < pre_win) = [];
locs_tri(locs_tri > size(EMG,1) - pre_win) = [];

[~, locs_bi] = findpeaks(EMG_smo(:,Bi_ch),'MinPeakHeight',0.2,'MinPeakDistance',80,'MinPeakProminence',0.2);
locs_bi(locs_bi < pre_win) = [];
locs_bi(locs_bi > size(EMG,1) - pre_win) = [];

num_events = length(locs_tri);
events_tri = zeros(num_events,pre_win*2+1);
events_bi = zeros(num_events,pre_win*2+1);
for i = 1:num_events
    events_tri(i,:) = EMG_smo(locs_tri(i)-pre_win:locs_tri(i)+pre_win,Tri_ch);
    events_bi(i,:) = EMG_smo(locs_tri(i)-pre_win:locs_tri(i)+pre_win,Bi_ch);
end
tri_peak_val = max(mean(events_tri(i,:)));

%Find modal duration from tri peak to next bi peak
durs_pos = zeros(num_events,1);
durs_neg = zeros(num_events,1);
for i = 1:num_events
    diffs = (locs_bi - locs_tri(i));
    diffs_pos = diffs(diffs>0);
    diffs_neg = diffs(diffs<0);
    if ~isempty(diffs_pos)
        durs_pos(i) = min(diffs_pos);
    end
    if ~isempty(diffs_neg)
        durs_neg(i) = max(diffs_neg);
    end
end
[n,bins_pos] = histc(durs_pos,50:5:pre_win); %peaks must be at least 50 ms later
[m,ind_pos] = max(n)
dur_pos_mode = 50+ind_pos*5-2.5
median_dur_pos = median(durs_pos)

[n,bins_neg] = histc(durs_neg,-pre_win:5:-50); %peaks must be at least 50 ms later
[m,ind_neg] = max(n)
dur_neg_mode = -pre_win+ind_neg*5-2.5
median_dur_neg = median(durs_neg)

modal_events = [];
width = 0;
while length(modal_events) < modal_event_set_size
    modal_pos_events = find(bins_pos <= ind_pos+width & bins_pos >= ind_pos-width);
    modal_neg_events = find(bins_neg <= ind_neg+width & bins_neg >= ind_neg-width);
    modal_events = intersect(modal_pos_events,modal_neg_events)
    width = width+1;
end
modal_tri_trace_mean = mean(events_tri(modal_events,:));
modal_tri_trace_sem = std(events_tri(modal_events,:))/sqrt(length(modal_events));

modal_bi_trace_mean = mean(events_bi(modal_events,:));
modal_bi_trace_sem = std(events_bi(modal_events,:))/sqrt(length(modal_events));

%Plot modal averages
figure; 
time = -pre_win:pre_win;
boundedline(time,modal_tri_trace_mean,modal_tri_trace_sem,'-b','transparency', 0.5)
hold on; 
boundedline(time,modal_bi_trace_mean,modal_bi_trace_sem,'-r','transparency', 0.5)
xlabel('Time from triceps peak (ms)')
ylabel('Normalized EMG magnitude (ms)')

%Load Ca data
load(CNMF_filename);
num_cells = size(A_or)
%Get row in which each cell was imaged


