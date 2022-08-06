%%%% Process CSV files with EMG, etc. and output from Simons processing
%%%According to Bruker, collection begins with first galvo rotation and
%%%does not skip any, such that only galvo rotations on the end of the time
%%%series should be ignored. Verify with Pockels blanking signal in the
%%%future
close all
clear
runtime = tic;

animal = 'D08';
data_path = 'C:\Users\jam5064\Desktop\Current Data\D08_imaging_data\TSeries\TSeries-11032017';
csv_filenames = {'TSeries-11032017-D08-001_Cycle00001_VoltageRecording_001.csv','TSeries-11032017-D08-002_Cycle00001_VoltageRecording_001.csv'};
Ca_filename = 'TSeries-11032017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
mph = 0.3; %MinPeakHeight

names = strsplit(csv_filenames{1},'-');
voltages_filename = [animal '_' char(names{2}) '_voltages'];
segments_filename = [animal '_' char(names{2}) '_segments'];
names = strsplit(Ca_filename,'-');
names2 = strsplit(char(names{4}),'.');
align_filename = [animal '_' char(names{2}) '_align_v11_' char(names2{1})];

preloaded_v = 1; % if data were already loaded, saves time
preloaded_beh = 1;
preloaded_Ca_ts = 1;

%Set Params
frames_averaged = 4;
num_rows_orig = 256;
num_frames_ind = 18000; %FOR INDIVIDUAL FILES
num_muscles = 4;
Y_galvo_col = 2; %column in csv of Y_galvo signal
Pockels_col = 3; %column in csv of Pockels blanking signal 
encoder_col = 4; %column in csv of encoder signal
solenoid_col = 5; %column in csv of solenoid signal
EMG_cols = 6:9; %columns of EMG signals in csv; change if one or more channels is dead
samp_rate = 10000; %original sampling rate. Probably 10000
lp = 40; %low pass filter freq (Hz) for EMG filtering
smoFiltSTD = 10; %STD of Gaussian for EMG filtering, in ms
sub_rate = 10; %rate at which to downsample data - to go from 10K to 1k, use 10
Bi_ch = 1; %EMG header position where the biceps signal comes from
Tri_ch = 2; %EMG header position where the triceps signal comes from
EDC_ch = 3; %EMG header position where the EDC signal comes from
PL_ch = 4; %EMG header position where the PL signal comes from
alt_event_dur = 0.3; % duration of alternation event epochs, in sec
modal_event_alt_size = 30;
inact_dur = 0.3; %in s, duration of inactivity segments
inact_event_size = 30;
cocon_pre_dur = 0.10; %duration pre triceps peak for cocon
cocon_post_dur = 0.20; %duration post triceps peak for cocon
cocon_event_size = 30;
cocon_peak_int = 100; % in samples
J = 10; %Number of iterations

%Initialize
next_fig_num = 1;

%%%%%%%%%%Read data
Pockels_full = [];
Y_galvo_full = [];
Y_galvo_unsub_full = [];
encoder_full = [];
solenoid_full = [];
EMG_full = [];

if ~preloaded_v
    for i = 1:size(csv_filenames,2)
        i
        cd(data_path)
        csv_readtime = tic;
        disp('LOADING CSV')
        Mat = csvread(csv_filenames{i},1,0) ;
        toc(csv_readtime)

        Pockels_unsub = Mat(:,Pockels_col);
        Y_galvo_unsub = Mat(:,Y_galvo_col);

        %Find the starting sample for scanning
        a = Y_galvo_unsub > 0;
        b = [diff(a); 0];
        frame_samps = find(b==-1);
        [~,samp_start] = max(Y_galvo_unsub(1:frame_samps(1)));
        
        %Find the end of the last collected frame
        last_collected_frame_samp = frame_samps(num_frames_ind*frames_averaged);
        c = unique(diff(frame_samps));
        [d,ind] = min(Y_galvo_unsub(last_collected_frame_samp:last_collected_frame_samp+c(1)));
        samp_end = last_collected_frame_samp+ind-1;

        %sanity check - make sure this jives with when the Pockels blanks for good
        a = find(Pockels_unsub>1);
        b = a(end);
        last_frame_to_pockels_off_time_ms = b-samp_end %Print this out for sanity check
        %should be much smaller than frame time, 16.67 ms
        
        %Downsample and separate signals
        Pockels = downsample(Mat(samp_start:samp_end,Pockels_col),sub_rate);
        Y_galvo = downsample(Mat(samp_start:samp_end,Y_galvo_col),sub_rate);
        Y_galvo_unsub = Mat(samp_start:samp_end,Y_galvo_col);
        encoder = downsample(Mat(samp_start:samp_end,encoder_col),sub_rate);
        solenoid = downsample(Mat(samp_start:samp_end,solenoid_col),sub_rate);
        EMG = downsample(Mat(samp_start:samp_end,EMG_cols),sub_rate);
        
        Pockels_full = [Pockels_full; Pockels];
        Y_galvo_full = [Y_galvo_full; Y_galvo];
        Y_galvo_unsub_full = [Y_galvo_unsub_full; Y_galvo_unsub];
        encoder_full = [encoder_full; encoder];
        solenoid_full = [solenoid_full; solenoid];
        EMG_full = [EMG_full; EMG];
    end
    
    %save downsampled csv stuff for easier reloading
    save(voltages_filename,'Pockels_full','Y_galvo_full','Y_galvo_unsub_full','encoder_full','solenoid_full','EMG_full')

else
    cd(data_path)
    load(voltages_filename);
end

Pockels = Pockels_full;
Y_galvo = Y_galvo_full;
Y_galvo_unsub = Y_galvo_unsub_full;
encoder = encoder_full; 
solenoid = solenoid_full;
EMG = EMG_full;
num_v_samps = size(EMG,1);

clear Pockels_full
clear Y_galvo_full
clear Y_galvo_unsub_full
clear encoder_full
clear solenoid_full
clear EMG_full
        
if ~preloaded_beh
    %Break the experiment up into epochs
    figure(next_fig_num); next_fig_num = next_fig_num+1; 
    a = smooth(diff(encoder),10,'sgolay',3);
    plot(a); hold on
    b = a > 0.01;
    wheel_crossings = find(diff(b)==-1);
    c = diff(wheel_crossings);
    [m,ind] = max(c);
    alt1_end = wheel_crossings(ind) + samp_rate/sub_rate; %add a one second buffer after the last wheel turn
    cocon_end = alt1_end + m - 2*samp_rate/sub_rate; %add a one second buffer before the first wheel turn
    plot([alt1_end cocon_end],[-2 -2],'ro')

    %%%%%%%%%%Process the EMG signals
    EMG_smo = 0*EMG;
    for i = 1:size(EMG,2)
        EMG_smo(:,i) = filterEMG(EMG(:,i),lp,smoFiltSTD);
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
    %    figure(next_fig_num); next_fig_num = next_fig_num+1; plot(EMG_smo(:,i))
    end

    %%%%%%%%%%Find events in alt session 1 (alt1)
    %May not want to find events in last 400 ms, according to apparent extra
    %voltage ts collection
    pre_win_alt = floor(alt_event_dur*samp_rate/sub_rate/2);
    [~, locs_tri1] = findpeaks(EMG_smo(:,Tri_ch),'MinPeakHeight',mph,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_tri1(locs_tri1 < pre_win_alt) = [];
    locs_tri1(locs_tri1 > alt1_end - pre_win_alt) = [];

    [~, locs_bi1] = findpeaks(EMG_smo(:,Bi_ch),'MinPeakHeight',mph,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_bi1(locs_bi1 < pre_win_alt) = [];
    locs_bi1(locs_bi1 > alt1_end - pre_win_alt) = [];

    %Find modal duration from tri peak to next bi peak
    durs_pos = zeros(length(locs_tri1),1);
    durs_neg = zeros(length(locs_tri1),1);
    for i = 1:length(locs_tri1)
        diffs = (locs_bi1 - locs_tri1(i));
        diffs_pos = diffs(diffs>0);
        diffs_neg = diffs(diffs<0);
        if ~isempty(diffs_pos)
            durs_pos(i) = min(diffs_pos);
        end
        if ~isempty(diffs_neg)
            durs_neg(i) = max(diffs_neg);
        end
    end
    [n,bins_pos] = histc(durs_pos,50:2.5:pre_win_alt); %peaks must be at least 50 ms later
    [m,ind_pos] = max(n);
    dur_pos_mode = 50+ind_pos*2.5-1.25;
    median_dur_pos = median(durs_pos);

    [n,bins_neg] = histc(durs_neg,-pre_win_alt:2.5:-50); %peaks must be at least 50 ms later
    [m,ind_neg] = max(n);
    dur_neg_mode = -pre_win_alt+ind_neg*2.5-1.25;
    median_dur_neg = median(durs_neg);

    modal_events_alt1 = [];
    width = 0;
    while length(modal_events_alt1) < modal_event_alt_size
        modal_pos_events = find(bins_pos <= ind_pos+width & bins_pos >= ind_pos-width);
        modal_neg_events = find(bins_neg <= ind_neg+width & bins_neg >= ind_neg-width);
        modal_events_alt1 = intersect(modal_pos_events,modal_neg_events);
        width = width+1;
    end

    num_alt1_events = length(modal_events_alt1);
    events_alt1_emg = zeros(num_alt1_events,pre_win_alt*2+1,4);
    for i = 1:num_alt1_events
        events_alt1_emg(i,:,1) = EMG_smo(locs_tri1(modal_events_alt1(i))-pre_win_alt:locs_tri1(modal_events_alt1(i))+pre_win_alt,Bi_ch);
        events_alt1_emg(i,:,2) = EMG_smo(locs_tri1(modal_events_alt1(i))-pre_win_alt:locs_tri1(modal_events_alt1(i))+pre_win_alt,Tri_ch);
        events_alt1_emg(i,:,3) = EMG_smo(locs_tri1(modal_events_alt1(i))-pre_win_alt:locs_tri1(modal_events_alt1(i))+pre_win_alt,EDC_ch);
        events_alt1_emg(i,:,4) = EMG_smo(locs_tri1(modal_events_alt1(i))-pre_win_alt:locs_tri1(modal_events_alt1(i))+pre_win_alt,PL_ch);
    end

    %Plot modal alt averages
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    time_alt = -pre_win_alt:pre_win_alt;
    boundedline(time_alt,mean(events_alt1_emg(:,:,2)),std(events_alt1_emg(:,:,2))/sqrt(num_alt1_events),'-b','transparency', 0.5)
    hold on; 
    boundedline(time_alt,mean(events_alt1_emg(:,:,1)),std(events_alt1_emg(:,:,1))/sqrt(num_alt1_events),'-r','transparency', 0.5)
    xlabel('Time from triceps peak (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Alternation 1 EMG')
    ylim([-0.05 1])
    hold off

    %%%%%%%%%%Find events in alt session 2 (alt2)
    %May not want to find events in last 400 ms, according to apparent extra
    %voltage ts collection
    pre_win_alt = floor(alt_event_dur*samp_rate/sub_rate/2);
    [~, locs_tri2] = findpeaks(EMG_smo(:,Tri_ch),'MinPeakHeight',mph,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_tri2(locs_tri2 < cocon_end) = [];
    locs_tri2(locs_tri2 > num_v_samps - pre_win_alt) = [];

    [~, locs_bi2] = findpeaks(EMG_smo(:,Bi_ch),'MinPeakHeight',mph,'MinPeakDistance',80,'MinPeakProminence',0.2);
    locs_bi2(locs_bi2 < cocon_end) = [];
    locs_bi2(locs_bi2 > num_v_samps - pre_win_alt) = [];

    %for each tri peak, check to see if there are bi peaks near by in both dirs
    elim = [];
    for i = 1:length(locs_tri1)
        diffs = (locs_bi1 - locs_tri1(i));
        diffs_pos = diffs(diffs>0);
        diffs_neg = diffs(diffs<0);
        dur_pos = Inf;
        if ~isempty(diffs_pos)
            dur_pos = min(diffs_pos);
        end
        dur_neg = Inf;
        if ~isempty(diffs_neg)
            dur_neg = max(diffs_neg);
        end
        if dur_pos > pre_win_alt || dur_neg < -pre_win_alt
            elim = [elim; i];
        end
    end
    locs_tri1(elim) = [];

    elim = [];
    for i = 1:length(locs_tri2)
        diffs = (locs_bi2 - locs_tri2(i));
        diffs_pos = diffs(diffs>0);
        diffs_neg = diffs(diffs<0);
        dur_pos = Inf;
        if ~isempty(diffs_pos)
            dur_pos = min(diffs_pos);
        end
        dur_neg = Inf;
        if ~isempty(diffs_neg)
            dur_neg = max(diffs_neg);
        end
        if dur_pos > pre_win_alt || dur_neg < -pre_win_alt
            elim = [elim; i];
        end
    end
    locs_tri2(elim) = [];

    num_alt2_events = num_alt1_events;

    %Now iterate back and forth
    mean_dists = zeros(2,J);
    for j = 1:J
        dists = zeros(length(locs_tri2),1);
        for i = 1:length(locs_tri2)
            d1 = pdist2(mean(events_alt1_emg(:,:,1)),EMG_smo(locs_tri2(i)-pre_win_alt:locs_tri2(i)+pre_win_alt,Bi_ch)');
            d2 = pdist2(mean(events_alt1_emg(:,:,2)),EMG_smo(locs_tri2(i)-pre_win_alt:locs_tri2(i)+pre_win_alt,Tri_ch)');
            dists(i) = sqrt(d1^2 + d2^2);
        end
        [~,lowest_dists_ind2] = sort(dists);
        mean_dists(1,j) = mean(dists(lowest_dists_ind2(1:num_alt2_events)));

        events_alt2_emg = zeros(num_alt2_events,pre_win_alt*2+1,4);
        for i = 1:num_alt2_events
            events_alt2_emg(i,:,1) = EMG_smo(locs_tri2(lowest_dists_ind2(i))-pre_win_alt:locs_tri2(lowest_dists_ind2(i))+pre_win_alt,Bi_ch);
            events_alt2_emg(i,:,2) = EMG_smo(locs_tri2(lowest_dists_ind2(i))-pre_win_alt:locs_tri2(lowest_dists_ind2(i))+pre_win_alt,Tri_ch);
            events_alt2_emg(i,:,3) = EMG_smo(locs_tri2(lowest_dists_ind2(i))-pre_win_alt:locs_tri2(lowest_dists_ind2(i))+pre_win_alt,EDC_ch);
            events_alt2_emg(i,:,4) = EMG_smo(locs_tri2(lowest_dists_ind2(i))-pre_win_alt:locs_tri2(lowest_dists_ind2(i))+pre_win_alt,PL_ch);
        end

        dists = zeros(length(locs_tri1),1);
        for i = 1:length(locs_tri1)
            d1 = pdist2(mean(events_alt2_emg(:,:,1)),EMG_smo(locs_tri1(i)-pre_win_alt:locs_tri1(i)+pre_win_alt,Bi_ch)');
            d2 = pdist2(mean(events_alt2_emg(:,:,2)),EMG_smo(locs_tri1(i)-pre_win_alt:locs_tri1(i)+pre_win_alt,Tri_ch)');
            dists(i) = sqrt(d1^2 + d2^2);
        end
        [~,lowest_dists_ind1] = sort(dists);
        mean_dists(2,j) = mean(dists(lowest_dists_ind1(1:num_alt1_events)));

        for i = 1:num_alt1_events
            events_alt1_emg(i,:,1) = EMG_smo(locs_tri1(lowest_dists_ind1(i))-pre_win_alt:locs_tri1(lowest_dists_ind1(i))+pre_win_alt,Bi_ch);
            events_alt1_emg(i,:,2) = EMG_smo(locs_tri1(lowest_dists_ind1(i))-pre_win_alt:locs_tri1(lowest_dists_ind1(i))+pre_win_alt,Tri_ch);
            events_alt1_emg(i,:,3) = EMG_smo(locs_tri1(lowest_dists_ind1(i))-pre_win_alt:locs_tri1(lowest_dists_ind1(i))+pre_win_alt,EDC_ch);
            events_alt1_emg(i,:,4) = EMG_smo(locs_tri1(lowest_dists_ind1(i))-pre_win_alt:locs_tri1(lowest_dists_ind1(i))+pre_win_alt,PL_ch);
        end

    end   

    figure(next_fig_num); next_fig_num = next_fig_num+1;
    plot(mean_dists(1,:),'r')
    hold on
    plot(mean_dists(2,:),'k')
    ylabel('Mean dists')

    %Plot modal alt averages
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    time_alt = -pre_win_alt:pre_win_alt;
    boundedline(time_alt,mean(events_alt1_emg(:,:,2)),std(events_alt1_emg(:,:,2))/sqrt(num_alt1_events),'-b','transparency', 0.5)
    hold on; 
    boundedline(time_alt,mean(events_alt1_emg(:,:,1)),std(events_alt1_emg(:,:,1))/sqrt(num_alt1_events),'-r','transparency', 0.5)
    xlabel('Time from triceps peak (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Alternation 1 EMG')
    ylim([-0.05 1])
    hold off

    %Plot modal alt averages
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    time_alt = -pre_win_alt:pre_win_alt;
    boundedline(time_alt,mean(events_alt2_emg(:,:,2)),std(events_alt2_emg(:,:,2))/sqrt(num_alt2_events),'-b','transparency', 0.5)
    hold on; 
    boundedline(time_alt,mean(events_alt2_emg(:,:,1)),std(events_alt2_emg(:,:,1))/sqrt(num_alt2_events),'-r','transparency', 0.5)
    xlabel('Time from triceps peak (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Alternation 2 EMG')
    ylim([-0.05 1])
    hold off

    %%%%%%%%%%Find Co-contraction events
    cocon_pre_samps = cocon_pre_dur*samp_rate/sub_rate;
    cocon_post_samps = cocon_post_dur*samp_rate/sub_rate;
    max_alt_peak_bi = max(mean(events_alt2_emg(:,1:end/2,1))); %onlt consider first biceps peak
    max_alt_peak_tri = max(mean(events_alt2_emg(:,:,2)));
    max_alt_peak_mean = mean([max_alt_peak_tri max_alt_peak_bi]);
    events_cocon_emg = zeros(cocon_event_size,cocon_pre_samps+cocon_post_samps+1,4);

    %Find all contiguous stretchs where both bi and tri are less than 0.1, find modal time to peak in tri,...
    a = EMG_smo(alt1_end:cocon_end,Tri_ch) < 0.1 & EMG_smo(alt1_end:cocon_end,Bi_ch) < 0.1;
    b = [diff(a); 0];
    drops = find(b==1);
    rises = find(b==-1);
    next_rise = zeros(length(drops)-1,1);
    for i = 1:length(drops)-1
        good_rises = find(rises > drops(i));
        next_rise(i) = rises(good_rises(1));
    end
    dr_pairs = [drops(1:end-1) next_rise];
    widths = diff(dr_pairs,1,2);

%     while mean([max(mean(events_cocon_emg(:,:,1))) max(mean(events_cocon_emg(:,:,2)))]) < max_alt_peak_mean-0.05 %just needs to be close
        crossings = next_rise(widths >= cocon_pre_samps);

        %Now remove those with peaks not high enough
        elim = [];
        for i = 1:length(crossings)
            if max(EMG_smo(alt1_end-1+crossings(i)+1:alt1_end-1+crossings(i)+cocon_peak_int,Bi_ch)) < mph || max(EMG_smo(alt1_end-1+crossings(i)+1:alt1_end-1+crossings(i)+cocon_peak_int,Tri_ch)) < mph
                elim = [elim; i];
            end
        end
        crossings(elim) = [];

        %Now find those with smallest pdist
        pdists = zeros(length(crossings),1);
        for i = 1:length(crossings)
            pdists(i) = pdist2(EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,Bi_ch)',EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,Tri_ch)');
        end
        [~,lowest_dists_ind_cocon] = sort(pdists);
        crossings = crossings(lowest_dists_ind_cocon(1:cocon_event_size));

        for i = 1:cocon_event_size
            events_cocon_emg(i,:,1) = EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,Bi_ch);
            events_cocon_emg(i,:,2) = EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,Tri_ch);
            events_cocon_emg(i,:,3) = EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,EDC_ch);
            events_cocon_emg(i,:,4) = EMG_smo(alt1_end-1+crossings(i)-cocon_pre_samps:alt1_end-1+crossings(i)+cocon_post_samps,PL_ch);
        end
        
%         mph = mph+0.01;

%     end

    %Plot cocon averages
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    time_cocon = -cocon_pre_samps:cocon_post_samps;
    boundedline(time_cocon,mean(events_cocon_emg(:,:,2)),std(events_cocon_emg(:,:,2))/sqrt(cocon_event_size),'-b','transparency', 0.5)
    hold on; 
    boundedline(time_cocon,mean(events_cocon_emg(:,:,1)),std(events_cocon_emg(:,:,1))/sqrt(cocon_event_size),'-r','transparency', 0.5)
    xlabel('Time (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Cocontraction EMG')
    ylim([-0.05 1])
    hold off

    %compute tial selection stats
    std_tri_A1 = mean(std(events_alt1_emg(:,:,2)))
    std_tri_A2 = mean(std(events_alt2_emg(:,:,2)))
    std_tri_CC = mean(std(events_cocon_emg(:,:,2)))

    std_bi_A1 = mean(std(events_alt1_emg(:,:,1)))
    std_bi_A2 = mean(std(events_alt2_emg(:,:,1)))
    std_bi_CC = mean(std(events_cocon_emg(:,:,1)))

    rmse_tri_A1A2 = sqrt(mean((mean(events_alt1_emg(:,:,2)) - mean(events_alt2_emg(:,:,2))).^2))
    rmse_bi_A1A2 = sqrt(mean((mean(events_alt1_emg(:,:,1)) - mean(events_alt2_emg(:,:,1))).^2))
    rmse_tribi_CC = sqrt(mean((mean(events_cocon_emg(:,:,1)) - mean(events_cocon_emg(:,:,2))).^2))

    peak_amp_std = std([max(mean(events_alt1_emg(:,:,1))) max(mean(events_alt1_emg(:,:,2))) max(mean(events_alt2_emg(:,:,1))) max(mean(events_alt2_emg(:,:,2))) max(mean(events_cocon_emg(:,:,1))) max(mean(events_cocon_emg(:,:,2)))])

    %%%%%%%%%%Find inactivity events
    inact_dur_samps = inact_dur*samp_rate/sub_rate;
    summed_EMG = [];
    wb=waitbar(0,'Finding inactive segments');
    for t = 1:10:size(EMG_smo,1)-inact_dur_samps
        summed_EMG = [summed_EMG; t sum(sum(EMG_smo(t:t+inact_dur_samps-1,:)))];
        waitbar(t/size(EMG_smo,1),wb,'Finding inactive segments');
    end
    close(wb)
    summed_EMG_sorted = sortrows(summed_EMG,2);

    %Now eliminate the overlapping chunks
    viols = 1;
    while ~isempty(viols)
        dummy = abs(bsxfun(@minus,summed_EMG_sorted(1:inact_event_size,1),summed_EMG_sorted(1:inact_event_size,1)')) < inact_dur_samps;
        [row,col] = find(triu(dummy,1));
        viols = unique(col);
        summed_EMG_sorted(viols,:) = [];
    end

    events_inact_emg = zeros(inact_event_size,inact_dur_samps,4);
    for i = 1:inact_event_size
        events_inact_emg(i,:,1) = EMG_smo(summed_EMG_sorted(i,1):summed_EMG_sorted(i,1)+inact_dur_samps-1,Bi_ch);         
        events_inact_emg(i,:,2) = EMG_smo(summed_EMG_sorted(i,1):summed_EMG_sorted(i,1)+inact_dur_samps-1,Tri_ch);
        events_inact_emg(i,:,3) = EMG_smo(summed_EMG_sorted(i,1):summed_EMG_sorted(i,1)+inact_dur_samps-1,EDC_ch);         
        events_inact_emg(i,:,4) = EMG_smo(summed_EMG_sorted(i,1):summed_EMG_sorted(i,1)+inact_dur_samps-1,PL_ch);
    end
    num_inact_events = inact_event_size;

    %Plot inact averages
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    time_inact = 1:inact_dur_samps;
    boundedline(time_inact,mean(events_inact_emg(:,:,2)),std(events_inact_emg(:,:,2))/sqrt(num_inact_events),'-b','transparency', 0.5)
    hold on; 
    boundedline(time_inact,mean(events_inact_emg(:,:,1)),std(events_inact_emg(:,:,1))/sqrt(num_inact_events),'-r','transparency', 0.5)
    boundedline(time_inact,mean(events_inact_emg(:,:,3)),std(events_inact_emg(:,:,3))/sqrt(num_inact_events),'-m','transparency', 0.5)
    boundedline(time_inact,mean(events_inact_emg(:,:,4)),std(events_inact_emg(:,:,4))/sqrt(num_inact_events),'-c','transparency', 0.5)
    xlabel('Time (ms)')
    ylabel('Normalized EMG magnitude (ms)')
    title('Inactivation EMG')
    ylim([-0.05 1])
    hold off

    toc(runtime)

    save(segments_filename,'EMG_smo','time_inact','events_inact_emg','time_cocon','events_cocon_emg','time_alt','events_alt1_emg','events_alt2_emg',...
            'locs_tri1','locs_tri2','lowest_dists_ind1','lowest_dists_ind2','alt1_end','crossings','summed_EMG_sorted','pre_win_alt',...
            'cocon_pre_samps','cocon_post_samps','inact_dur_samps');
else
    load(segments_filename)
end

%%%%%%%%%%Now make time series for the cells
%Load Ca data
load(Ca_filename);
num_cells = size(A,2);

if ~preloaded_Ca_ts
    num_rows = 256-length(rows);
    num_cols = 256-length(cols);
    num_frames = size(R,2);

    %Assuming that first galvo rotation in ts is first collected frame, and no galvo rotations are skipped during collection NEED TO VERIFY - need xml data
    %also using "time" here as a proxy for samples - this is ok as long as the 1 sample = 1 ms convention is maintained
    peak_vals = findpeaks(Y_galvo_unsub);
    min_peak = min([Y_galvo_unsub(1); peak_vals]);
    trough_vals = -findpeaks(-Y_galvo_unsub);
    max_trough = max([Y_galvo_unsub(end); trough_vals]);
    galvo_bnds = quantile(max_trough:0.001:min_peak,[1/(2*num_rows) 1-1/(2*num_rows)]); %avoid the ends where Y_galvo may hit extremes of range not seen during scanning at higher zone

    Cell_ts = zeros(size(Y_galvo,1),num_cells,3);

    wb=waitbar(0,'Interpolating cellular time series');
    for i = 1:num_cells
        Mat = reshape(A(:,i),num_rows,num_cols,1); 
        cOM = centerOfMass(Mat);
        cell_row_pos = cOM(1);
        row_galvo_pos = galvo_bnds(2)- abs(diff(galvo_bnds)/(num_rows-1))*(cell_row_pos-1);
        a = Y_galvo_unsub > row_galvo_pos;
        b = diff(a);
        row_times_full = find(b==-1);
        row_times = round(row_times_full/sub_rate); %Get back to downsampled time points
        if length(row_times)~= frames_averaged*num_frames
            disp('ISSUE: length(row_times)~= frames_averaged*num_frames')
        end
        c = diff(row_times);
        d = unique(c);
        half_frame_time = floor(d(1)/2);
        for j = 1:num_frames
            bnds = [max([1 row_times(frames_averaged*(j-1)+1)-half_frame_time])  min([size(Y_galvo,1) row_times(frames_averaged*j)+half_frame_time])];
            Cell_ts(bnds(1):bnds(2),i,1) = Rsub(i,j);
            Cell_ts(bnds(1):bnds(2),i,2) = Rdcv(i,j); 
            Cell_ts(bnds(1):bnds(2),i,3) = S(i,j)*frames_averaged; %I think this is right - sp is computed as spikes per time bin/frame?
        end
        %Smooth Cell time series to be equivalent to EMG, getting rid of rough edges
        for j = 1:3
            Cell_ts(:,i,j) = filterGauss(Cell_ts(:,i,j),smoFiltSTD);
            if j==2
                Cell_ts(:,i,j) = zscore(Cell_ts(:,i,j)); %Z-score the Ca signal so PCs are mag independent
            end
        end

        waitbar(i/num_cells,wb,'Interpolating cellular time series');
    end
    
    %%%Compute Cell Corrs with muscles
    cell_corrs = corr(Cell_ts(:,:,1),EMG_smo);
    cell_corrs(:,:,2) = corr(Cell_ts(:,:,2),EMG_smo); 
    cell_corrs(:,:,3) = corr(Cell_ts(:,:,3),EMG_smo);

    num_alt1_events = size(events_alt1_emg,1);
    events_alt1_cell = zeros(num_alt1_events,pre_win_alt*2+1,num_cells,3);

    num_alt2_events = size(events_alt2_emg,1);
    events_alt2_cell = zeros(num_alt2_events,pre_win_alt*2+1,num_cells,3);

    num_cocon_events = size(events_cocon_emg,1);
    events_cocon_cell = zeros(num_cocon_events,cocon_pre_samps+cocon_post_samps+1,num_cells,3);

    num_inact_events = size(events_inact_emg,1);
    events_inact_cell = zeros(num_inact_events,inact_dur_samps,num_cells,3);

    for i = 1:num_cells
        for k = 1:3
            for j = 1:num_alt1_events
                events_alt1_cell(j,:,i,k) = Cell_ts(locs_tri1(lowest_dists_ind1(j))-pre_win_alt:locs_tri1(lowest_dists_ind1(j))+pre_win_alt,i,k);
            end
            for j = 1:num_alt2_events
                events_alt2_cell(j,:,i,k) = Cell_ts(locs_tri2(lowest_dists_ind2(j))-pre_win_alt:locs_tri2(lowest_dists_ind2(j))+pre_win_alt,i,k);
            end
            for j = 1:num_cocon_events
                events_cocon_cell(j,:,i,k) = Cell_ts(alt1_end-1+crossings(j)-cocon_pre_samps:alt1_end-1+crossings(j)+cocon_post_samps,i,k);
            end
            for j = 1:num_inact_events
                events_inact_cell(j,:,i,k) = Cell_ts(summed_EMG_sorted(j,1):summed_EMG_sorted(j,1)+inact_dur_samps-1,i,k);
            end
        end
        waitbar(i/num_cells,wb,'Collecting cellular time series segments');
    end
    close(wb)

    save(Ca_filename,'cell_corrs','events_alt1_cell','events_alt2_cell','events_cocon_cell','events_inact_cell','-append')
end

alt1_means = zeros(pre_win_alt*2+1,num_cells,3);
alt2_means = zeros(pre_win_alt*2+1,num_cells,3);
cocon_means = zeros(cocon_pre_samps+cocon_post_samps+1,num_cells,3);
inact_means = zeros(inact_dur_samps,num_cells,3);

for i = 1:num_cells
    for j = 1:3
        alt1_means(:,i,j) = mean(events_alt1_cell(:,:,i,j));
        alt2_means(:,i,j) = mean(events_alt2_cell(:,:,i,j));
        cocon_means(:,i,j) = mean(events_cocon_cell(:,:,i,j));
        inact_means(:,i,j) = mean(events_inact_cell(:,:,i,j));
    end
end

%One for method for finding cells to keep  - keep cells that make
%significant contributions to top PCs
num_stds = 6;
noise_comp_cuttoff = round(num_cells/2);
align_pc_comp = zeros(2,3,4,3);
good_cells_full = [];
var_exp = cell(2,4,3);
for PCAdim = [5 6 7]
    %calcium - good cells define with spiking tho
    coeff_sp_alt1 = pca(alt1_means(:,:,3));
    noise_comps = abs(coeff_sp_alt1(:,end-noise_comp_cuttoff:end));
    cell_comps_mean = mean(noise_comps,2);
    cell_comps_std = std(noise_comps,0,2);
    comp_alt = repmat(cell_comps_mean + num_stds*cell_comps_std,1,PCAdim);
    
    coeff_sp_cocon = pca(cocon_means(:,:,3));
    noise_comps = abs(coeff_sp_cocon(:,end-noise_comp_cuttoff:end));
    cell_comps_mean = mean(noise_comps,2);
    cell_comps_std = std(noise_comps,0,2);
    comp_cocon = repmat(cell_comps_mean + num_stds*cell_comps_std,1,PCAdim);

    good_cells = find(sum(abs(coeff_sp_alt1(:,PCAdim))>comp_alt,2) + sum(abs(coeff_sp_cocon(:,PCAdim))>comp_cocon,2));
    good_cells_full(PCAdim-4).ca = good_cells;
    
    [coeff_ca_alt1,~,~,~,var_exp{1,1,PCAdim-4}] = pca(alt1_means(:,good_cells,2));
    [coeff_ca_cocon,~,~,~,var_exp{1,2,PCAdim-4}] = pca(cocon_means(:,good_cells,2));
    [coeff_ca_alt2,~,~,~,var_exp{1,3,PCAdim-4}] = pca(alt2_means(:,good_cells,2));
    [coeff_ca_inact,~,~,~,var_exp{1,4,PCAdim-4}] = pca(inact_means(:,good_cells,2));
    
    align_alt1alt2 = (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt2_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim))) / (trace(coeff_ca_alt2(:,1:PCAdim)'*cov(alt2_means(:,good_cells,2))*coeff_ca_alt2(:,1:PCAdim)));
    align_alt2alt1 = (trace(coeff_ca_alt2(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt2(:,1:PCAdim))) / (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim)));
    align_pc_comp(1,3,1,PCAdim-4) = mean([align_alt1alt2 align_alt2alt1]);
    
    align_alt1cocon = (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(cocon_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim))) / (trace(coeff_ca_cocon(:,1:PCAdim)'*cov(cocon_means(:,good_cells,2))*coeff_ca_cocon(:,1:PCAdim)));
    align_coconalt1 = (trace(coeff_ca_cocon(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_cocon(:,1:PCAdim))) / (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim)));
    align_pc_comp(1,3,2,PCAdim-4) = mean([align_alt1cocon align_coconalt1]);
    
    align_alt1inact = (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(inact_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim))) / (trace(coeff_ca_inact(:,1:PCAdim)'*cov(inact_means(:,good_cells,2))*coeff_ca_inact(:,1:PCAdim)));
    align_inactalt1 = (trace(coeff_ca_inact(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_inact(:,1:PCAdim))) / (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim)));
    align_pc_comp(1,3,3,PCAdim-4) = mean([align_alt1inact align_inactalt1]);
        
    align_alt1alt1 = (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim))) / (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim)));
    align_alt1alt1b = (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim))) / (trace(coeff_ca_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,2))*coeff_ca_alt1(:,1:PCAdim)));
    align_pc_comp(1,3,4,PCAdim-4) = mean([align_alt1alt1 align_alt1alt1b]);
    
    good_cells_full(PCAdim-4).sp = good_cells;
    
    [coeff_sp_alt1,~,~,~,var_exp{2,1,PCAdim-4}] = pca(alt1_means(:,good_cells,3));
    [coeff_sp_cocon,~,~,~,var_exp{2,2,PCAdim-4}] = pca(cocon_means(:,good_cells,3));
    [coeff_sp_alt2,~,~,~,var_exp{2,3,PCAdim-4}] = pca(alt2_means(:,good_cells,3));
    [coeff_sp_inact,~,~,~,var_exp{2,4,PCAdim-4}] = pca(inact_means(:,good_cells,3));
    
    align_alt1alt2 = (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt2_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim))) / (trace(coeff_sp_alt2(:,1:PCAdim)'*cov(alt2_means(:,good_cells,3))*coeff_sp_alt2(:,1:PCAdim)));
    align_alt2alt1 = (trace(coeff_sp_alt2(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt2(:,1:PCAdim))) / (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim)));
    align_pc_comp(2,3,1,PCAdim-4) = mean([align_alt1alt2 align_alt2alt1]);
    
    align_alt1cocon = (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(cocon_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim))) / (trace(coeff_sp_cocon(:,1:PCAdim)'*cov(cocon_means(:,good_cells,3))*coeff_sp_cocon(:,1:PCAdim)));
    align_coconalt1 = (trace(coeff_sp_cocon(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_cocon(:,1:PCAdim))) / (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim)));
    align_pc_comp(2,3,2,PCAdim-4) = mean([align_alt1cocon align_coconalt1]);
    
    align_alt1inact = (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(inact_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim))) / (trace(coeff_sp_inact(:,1:PCAdim)'*cov(inact_means(:,good_cells,3))*coeff_sp_inact(:,1:PCAdim)));
    align_inactalt1 = (trace(coeff_sp_inact(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_inact(:,1:PCAdim))) / (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim)));
    align_pc_comp(2,3,3,PCAdim-4) = mean([align_alt1inact align_inactalt1]);
        
    align_alt1alt1 = (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim))) / (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim)));
    align_alt1alt1b = (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim))) / (trace(coeff_sp_alt1(:,1:PCAdim)'*cov(alt1_means(:,good_cells,3))*coeff_sp_alt1(:,1:PCAdim)));
    align_pc_comp(2,3,4,PCAdim-4) = mean([align_alt1alt1 align_alt1alt1b]);

end

save(align_filename,'align_pc_comp','good_cells_full','var_exp')
