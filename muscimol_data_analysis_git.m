% analysis of muscimol data (from the ns5 files) 

%% below is run and saved, open the files from
% 'F:\cohort5_alt_behavior_and_musc' and 'F:\cohort5_cc_behavior_and_musc'
% and run code from COCON LEARNING

%%
% %% define variables
% 
% dur=20; % cut all sessions to first 20 minutes
% cond_key_alt=struct; %A37-39 are same experimental group, A40-41 other group
% %A37
% cond_key_alt(14).A37= 'baseline'; cond_key_alt(15).A37= 'contra musc';
% cond_key_alt(16).A37= 'wash contra musc'; cond_key_alt(17).A37= 'contra sal';
% cond_key_alt(18).A37= 'wash contra sal'; cond_key_alt(19).A37= 'bilat musc';
% cond_key_alt(20).A37= 'wash bilat musc'; cond_key_alt(21).A37= 'bilat sal';
% cond_key_alt(22).A37= 'wash bilat sal';
% %A38
% cond_key_alt(14).A38= 'baseline'; cond_key_alt(15).A38= 'contra musc';
% cond_key_alt(16).A38= 'wash contra musc'; cond_key_alt(17).A38= 'contra sal';
% cond_key_alt(18).A38= 'wash contra sal'; cond_key_alt(19).A38= 'bilat musc';
% cond_key_alt(20).A38= 'wash bilat musc'; cond_key_alt(21).A38= 'bilat sal';
% cond_key_alt(22).A38= 'wash bilat sal';
% %A39
% cond_key_alt(14).A39= 'baseline'; cond_key_alt(15).A39= 'contra musc';
% cond_key_alt(16).A39= 'wash contra musc'; cond_key_alt(17).A39= 'contra sal';
% cond_key_alt(18).A39= 'wash contra sal'; cond_key_alt(19).A39= 'bilat musc';
% cond_key_alt(20).A39= 'wash bilat musc'; cond_key_alt(21).A39= 'bilat sal';
% cond_key_alt(22).A39= 'wash bilat sal';
% %A40
% cond_key_alt(14).A40= 'baseline'; cond_key_alt(15).A40= 'contra sal';
% cond_key_alt(16).A40= 'wash contra sal'; cond_key_alt(17).A40= 'contra musc';
% cond_key_alt(18).A40= 'wash contra musc'; cond_key_alt(19).A40= 'bilat sal';
% cond_key_alt(20).A40= 'wash bilat sal'; cond_key_alt(21).A40= 'bilat musc';
% cond_key_alt(22).A40= 'wash bilat musc';
% %A41
% cond_key_alt(14).A41= 'baseline'; cond_key_alt(15).A41= 'contra sal';
% cond_key_alt(16).A41= 'wash contra sal'; cond_key_alt(17).A41= 'contra musc';
% cond_key_alt(18).A41= 'wash contra musc'; cond_key_alt(19).A41= 'bilat sal';
% cond_key_alt(20).A41= 'wash bilat sal'; cond_key_alt(21).A41= 'bilat musc';
% cond_key_alt(22).A41= 'wash bilat musc';
% 
% cd('F:\cohort5_alt_behavior_and_musc');
% muscimol_alt_data=struct;
% mouse_IDs={'A37','A38','A39','A40','A41'};
% 
% for k=1:length(mouse_IDs)
%     mouse_IDs{k}
%     cd(strcat(mouse_IDs{k},'_alt_training'));
%     files=dir;
%     files(1:2)=[];
%     alt_data=struct;
%     
%     for ii=1:length(files)
%         ii
%         filename=files(ii).name;
%         session=ii;
%         openNSx(filename,'read');
%         
%         %decimate all data to 1kHz
%         spdel=decimate(double(NS5.Data(1,:)),30);
%         pec=decimate(double(NS5.Data(2,:)),30);
%         triM=decimate(double(NS5.Data(3,:)),30);
%         bi=decimate(double(NS5.Data(5,:)),30);
%         edc=decimate(double(NS5.Data(7,:)),30);
%         pal=decimate(double(NS5.Data(8,:)),30);
%         
%         enc=downsample(double(NS5.Data(9,:)),30);
%         %lds=downsample(double(NS5.Data(11,:)),30);
%         
%         % convert to mV
%         scale_factor=100*(2^15-1)/5000; % 100x to account for gain on amp. 
%         % 2^15-1 as the digital range mapped onto analog range (-5k to 5k mV).
%         spdel=spdel/scale_factor;
%         pec=pec/scale_factor;
%         triM=triM/scale_factor;
%         bi=bi/scale_factor;
%         edc=edc/scale_factor;
%         pal=pal/scale_factor;
%         
%         % filter emg
%         sd=10;
%         samp_rate=1000;
%         spdelSm = filterEMG_vSR(spdel, 40, sd, samp_rate);
%         pecSm = filterEMG_vSR(pec, 40, sd, samp_rate);
%         triMSm = filterEMG_vSR(triM, 40, sd, samp_rate);
%         biSm = filterEMG_vSR(bi, 40, sd, samp_rate);
%         edcSm = filterEMG_vSR(edc, 40, sd, samp_rate);
%         palSm = filterEMG_vSR(pal, 40, sd, samp_rate);
%         
%         % process encoder
%         thr_encdiff=0.4;
%         minEnc=min(enc);
%         maxEnc=max(enc);
%         enc=((enc-minEnc)/(maxEnc-minEnc))*3.599; %map to 0-3.599
%         enc=3.6-enc;
%         diffEnc=diff(enc);
%         jumps_neg_locs=find(diffEnc<-thr_encdiff);
%         jumps_neg_vals=diffEnc(jumps_neg_locs);
%         jumps_pos_locs=find(diffEnc>thr_encdiff);
%         jumps_pos_vals=diffEnc(jumps_pos_locs);
%         for i=1:length(jumps_neg_locs)
%             enc=[enc(1:jumps_neg_locs(i)) enc(jumps_neg_locs(i)+1:end)-jumps_neg_vals(i)];
%         end
%         for i=1:length(jumps_pos_locs)
%             enc=[enc(1:jumps_pos_locs(i)) enc(jumps_pos_locs(i)+1:end)-jumps_pos_vals(i)];
%         end
%         dist_cm=enc*(20/3.6); %wheel is 20 cm in perimeter
%         enc_vel=filterGauss(diff(dist_cm),20); %here 20 ms = 20 samples, size of filter
%         enc_vel=enc_vel*1000; % convert to cm/s
%         
%         % process solenoid reward
%         sol=(double(NS5.Data(10,:)));
%         sol(sol<30000)=0;
%         sol_diff=diff(sol);
%         sol_diff(sol_diff<10000)=0;
%         rew_vec=decimate_spiketrain(sol_diff,30); %vector of binary vals for rew
%         
%         % log data into struct
%         alt_data(session).original_dur_minutes=(length(dist_cm)/1000)/60; % in minutes
%         alt_data(session).minEnc=minEnc;
%         alt_data(session).maxEnc=maxEnc;
%         end_ind=20*60*1000; %20 minutes into ms
%         alt_data(session).rew_vec=rew_vec(1:end_ind);
%         alt_data(session).dist_cm=dist_cm(1:end_ind);
%         alt_data(session).enc_vel=enc_vel(1:end_ind);
%         alt_data(session).spdel=spdelSm(1:end_ind);
%         alt_data(session).pec=pecSm(1:end_ind);
%         alt_data(session).triM=triMSm(1:end_ind);
%         alt_data(session).bi=biSm(1:end_ind);
%         alt_data(session).edc=edcSm(1:end_ind);
%         alt_data(session).pal=palSm(1:end_ind);
%         alt_data(session).velMean20min=mean(enc_vel);
%         alt_data(session).velMax20min=max(enc_vel);
%         
%         emg_mat=[sort(alt_data(session).spdel); sort(alt_data(session).pec); sort(alt_data(session).triM); ...
%             sort(alt_data(session).bi); sort(alt_data(session).edc); sort(alt_data(session).pal)];
%         alt_data(session).emg_sorted_mean=mean(emg_mat);
%         alt_data(session).emg_sorted_sem=std(emg_mat)/sqrt(6);
%     end
%     
%     cd ..
%     muscimol_alt_data.(mouse_IDs{k})=alt_data;
%     
% end
% 
% 
% % get CCI
% for i=1:length(mouse_IDs)
%     i
%     % get normalization values from alt baseline sessions
%     mouse_num=mouse_IDs{i};
%     baseline_day1=12;
%     baseline_day2=13;
%     baseline_day3=14;
%     tri_percentile_low1=prctile(muscimol_alt_data.(mouse_num)(baseline_day1).triM,0.01);
%     tri_percentile_low2=prctile(muscimol_alt_data.(mouse_num)(baseline_day2).triM,0.01);
%     tri_percentile_low3=prctile(muscimol_alt_data.(mouse_num)(baseline_day3).triM,0.01);
%     tri_percentile_low=mean([tri_percentile_low1 tri_percentile_low2 tri_percentile_low3]);
%     tri_percentile_high1=prctile(muscimol_alt_data.(mouse_num)(baseline_day1).triM,99.99);
%     tri_percentile_high2=prctile(muscimol_alt_data.(mouse_num)(baseline_day2).triM,99.99);
%     tri_percentile_high3=prctile(muscimol_alt_data.(mouse_num)(baseline_day3).triM,99.99);
%     tri_percentile_high=mean([tri_percentile_high1 tri_percentile_high2 tri_percentile_high3]);
%     bi_percentile_low1=prctile(muscimol_alt_data.(mouse_num)(baseline_day1).bi,0.01);
%     bi_percentile_low2=prctile(muscimol_alt_data.(mouse_num)(baseline_day2).bi,0.01);
%     bi_percentile_low3=prctile(muscimol_alt_data.(mouse_num)(baseline_day3).bi,0.01);
%     bi_percentile_low=mean([bi_percentile_low1 bi_percentile_low2  bi_percentile_low3]);
%     bi_percentile_high1=prctile(muscimol_alt_data.(mouse_num)(baseline_day1).bi,99.99);
%     bi_percentile_high2=prctile(muscimol_alt_data.(mouse_num)(baseline_day2).bi,99.99);
%     bi_percentile_high3=prctile(muscimol_alt_data.(mouse_num)(baseline_day3).bi,99.99);
%     bi_percentile_high=mean([bi_percentile_high1 bi_percentile_high2  bi_percentile_high3]);
% 
%     for session = 1:22
%         % normalize muscle traces
%         tri_norm=(muscimol_alt_data.(mouse_num)(session).triM-tri_percentile_low)/(tri_percentile_high-tri_percentile_low);
%         tri_norm(tri_norm<0)=0;
%         muscimol_alt_data.(mouse_num)(session).tri_norm=tri_norm;
%         bi_norm=(muscimol_alt_data.(mouse_num)(session).bi-bi_percentile_low)/(bi_percentile_high-bi_percentile_low);
%         bi_norm(bi_norm<0)=0;
%         muscimol_alt_data.(mouse_num)(session).bi_norm=bi_norm;
%         
%         % calculate CCI
%         CCI=nan*ones(1,length(muscimol_alt_data.(mouse_num)(session).bi_norm));
%         for k=1:length(CCI)
%             if muscimol_alt_data.(mouse_num)(session).bi_norm(k) > muscimol_alt_data.(mouse_num)(session).tri_norm(k)
%                 CCI(k)= (muscimol_alt_data.(mouse_num)(session).tri_norm(k)/muscimol_alt_data.(mouse_num)(session).bi_norm(k))*...
%                 (muscimol_alt_data.(mouse_num)(session).bi_norm(k)+muscimol_alt_data.(mouse_num)(session).tri_norm(k));
%             elseif muscimol_alt_data.(mouse_num)(session).tri_norm(k) > muscimol_alt_data.(mouse_num)(session).bi_norm(k)    
%             CCI(k)= (muscimol_alt_data.(mouse_num)(session).bi_norm(k)/muscimol_alt_data.(mouse_num)(session).tri_norm(k))*...
%                 (muscimol_alt_data.(mouse_num)(session).bi_norm(k)+muscimol_alt_data.(mouse_num)(session).tri_norm(k));
%             else %if they're the same, ratio will be 1
%                 CCI(k)= (muscimol_alt_data.(mouse_num)(session).bi_norm(k)/muscimol_alt_data.(mouse_num)(session).tri_norm(k))*...
%                 (muscimol_alt_data.(mouse_num)(session).bi_norm(k)+muscimol_alt_data.(mouse_num)(session).tri_norm(k));
%             end
%         end
%         muscimol_alt_data.(mouse_num)(session).CCI=CCI;
%         muscimol_alt_data.(mouse_num)(session).CCI_mean=mean(CCI);
%         muscimol_alt_data.(mouse_num)(session).CCI_median=median(CCI);
%         muscimol_alt_data.(mouse_num)(session).CCI_pmax=prctile(CCI,99.99); 
%     end
% end
% 
% %save data
% save('muscimol_alt_data','muscimol_alt_data','-v7.3');
% cd('F:');
% 
% %% cocon
% 
% cond_key_cc=struct; %A37-39 are same, A40-41 are same
% %A37
% cond_key_cc(15).A37= 'baseline'; cond_key_cc(16).A37= 'contra musc';
% cond_key_cc(17).A37= 'wash contra musc'; cond_key_cc(18).A37= 'contra sal';
% cond_key_cc(19).A37= 'wash contra sal'; cond_key_cc(20).A37= 'bilat musc';
% cond_key_cc(21).A37= 'wash bilat musc'; cond_key_cc(22).A37= 'bilat sal';
% cond_key_cc(23).A37= 'wash bilat sal'; cond_key_cc(24).A37= 'wash bilat sal redo';
% % A37's session 23 has shifted waterport placement, session was re-run (session 24)
% %A38
% cond_key_cc(15).A38= 'baseline'; cond_key_cc(16).A38= 'contra musc';
% cond_key_cc(17).A38= 'wash contra musc'; cond_key_cc(18).A38= 'contra sal';
% cond_key_cc(19).A38= 'wash contra sal'; cond_key_cc(20).A38= 'bilat musc';
% cond_key_cc(21).A38= 'wash bilat musc'; cond_key_cc(22).A38= 'bilat sal';
% cond_key_cc(23).A38= 'wash bilat sal';
% %A39
% cond_key_cc(15).A39= 'baseline'; cond_key_cc(16).A39= 'contra musc';
% cond_key_cc(17).A39= 'wash contra musc'; cond_key_cc(18).A39= 'contra sal';
% cond_key_cc(19).A39= 'wash contra sal'; cond_key_cc(20).A39= 'bilat musc';
% cond_key_cc(21).A39= 'wash bilat musc'; cond_key_cc(22).A39= 'bilat sal';
% cond_key_cc(23).A39= 'wash bilat sal';
% %A40
% cond_key_cc(15).A40= 'baseline'; cond_key_cc(16).A40= 'contra sal';
% cond_key_cc(17).A40= 'wash contra sal'; cond_key_cc(18).A40= 'contra musc';
% cond_key_cc(19).A40= 'wash contra musc'; cond_key_cc(20).A40= 'bilat sal';
% cond_key_cc(21).A40= 'wash bilat sal'; cond_key_cc(22).A40= 'bilat musc';
% cond_key_cc(23).A40= 'wash bilat musc';
% %A41
% cond_key_cc(15).A41= 'baseline'; cond_key_cc(16).A41= 'contra sal';
% cond_key_cc(17).A41= 'wash contra sal'; cond_key_cc(18).A41= 'contra musc';
% cond_key_cc(19).A41= 'wash contra musc'; cond_key_cc(20).A41= 'bilat sal';
% cond_key_cc(21).A41= 'wash bilat sal'; cond_key_cc(22).A41= 'bilat musc';
% cond_key_cc(23).A41= 'wash bilat musc';
% 
% cd('F:\cohort5_cocon_behavior_and_musc');
% muscimol_cc_data=struct;
% mouse_IDs={'A37','A38','A39','A40','A41'};
% 
% for k=1:length(mouse_IDs)
%     mouse_IDs{k}
%     cd(strcat(mouse_IDs{k},'_cc_training'));
%     files=dir;
%     files(1:2)=[];
%     
%      for ii=1:length(files)
%         ii
%         filename=files(ii).name;
%         session=ii;
%         openNSx(filename,'read');
%         
%         %downsample all data to 1kHz
%         spdel=decimate(double(NS5.Data(1,:)),30);
%         pec=decimate(double(NS5.Data(2,:)),30);
%         triM=decimate(double(NS5.Data(3,:)),30);
%         bi=decimate(double(NS5.Data(5,:)),30);
%         edc=decimate(double(NS5.Data(7,:)),30);
%         pal=decimate(double(NS5.Data(8,:)),30);
%         %enc=downsample(double(NS5.Data(9,:)),30);
%         %lds=downsample(double(NS5.Data(11,:)),30);
%         
%         % convert  to mV
%         scale_factor=100*(2^15-1)/5000; % 100x to account for gain on amp. 
%         % 2^15-1 as the digital range mapped onto analog range (-5k to 5k mV).
%         spdel=spdel/scale_factor;
%         pec=pec/scale_factor;
%         triM=triM/scale_factor;
%         bi=bi/scale_factor;
%         edc=edc/scale_factor;
%         pal=pal/scale_factor;
%         
%         % filter emg
%         sd=10;
%         samp_rate=1000;
%         spdelSm = filterEMG_vSR(spdel, 40, sd, samp_rate);
%         pecSm = filterEMG_vSR(pec, 40, sd, samp_rate);
%         triMSm = filterEMG_vSR(triM, 40, sd, samp_rate);
%         biSm = filterEMG_vSR(bi, 40, sd, samp_rate);
%         edcSm = filterEMG_vSR(edc, 40, sd, samp_rate);
%         palSm = filterEMG_vSR(pal, 40, sd, samp_rate);
%         
%         % process solenoid reward
%         sol=(double(NS5.Data(10,:)));
%         sol(sol<30000)=0;
%         sol_diff=diff(sol);
%         sol_diff(sol_diff<10000)=0;
%         rew_vec=decimate_spiketrain(sol_diff,30); %vector of binary vals for rew
%         
%         % log data into struct
%         cc_data(session).original_dur_minutes=(length(rew_vec)/1000)/60; % in minutes
%         end_ind=20*60*1000; %20 minutes into ms
%         cc_data(session).rew_vec=rew_vec(1:end_ind);
%         cc_data(session).spdel=spdelSm(1:end_ind);
%         cc_data(session).pec=pecSm(1:end_ind);
%         cc_data(session).triM=triMSm(1:end_ind);
%         cc_data(session).bi=biSm(1:end_ind);
%         cc_data(session).edc=edcSm(1:end_ind);
%         cc_data(session).pal=palSm(1:end_ind);
%         
%         corr_val=corrcoef(cc_data(session).triM,cc_data(session).bi);
%         cc_data(session).tot_tribi_corr = corr_val(1,2);
%         
%     end
%     
%     cd ..
%     muscimol_cc_data.(mouse_IDs{k})=cc_data;
% end
% 
% for i=1:length(mouse_IDs)
%     i
%     % get normalization values
%     mouse_num=mouse_IDs{i};
%     baseline_day1=13;
%     baseline_day2=14;
%     baseline_day3=15;
%     tri_percentile_low1=prctile(muscimol_cc_data.(mouse_num)(baseline_day1).triM,0.01);
%     tri_percentile_low2=prctile(muscimol_cc_data.(mouse_num)(baseline_day2).triM,0.01);
%     tri_percentile_low3=prctile(muscimol_cc_data.(mouse_num)(baseline_day3).triM,0.01);
%     tri_percentile_low=mean([tri_percentile_low1 tri_percentile_low2 tri_percentile_low3]);
%     tri_percentile_high1=prctile(muscimol_cc_data.(mouse_num)(baseline_day1).triM,99.99);
%     tri_percentile_high2=prctile(muscimol_cc_data.(mouse_num)(baseline_day2).triM,99.99);
%     tri_percentile_high3=prctile(muscimol_cc_data.(mouse_num)(baseline_day3).triM,99.99);
%     tri_percentile_high=mean([tri_percentile_high1 tri_percentile_high2 tri_percentile_high3]);
%     bi_percentile_low1=prctile(muscimol_cc_data.(mouse_num)(baseline_day1).bi,0.01);
%     bi_percentile_low2=prctile(muscimol_cc_data.(mouse_num)(baseline_day2).bi,0.01);
%     bi_percentile_low3=prctile(muscimol_cc_data.(mouse_num)(baseline_day3).bi,0.01);
%     bi_percentile_low=mean([bi_percentile_low1 bi_percentile_low2 bi_percentile_low3]);
%     bi_percentile_high1=prctile(muscimol_cc_data.(mouse_num)(baseline_day1).bi,99.99);
%     bi_percentile_high2=prctile(muscimol_cc_data.(mouse_num)(baseline_day2).bi,99.99);
%     bi_percentile_high3=prctile(muscimol_cc_data.(mouse_num)(baseline_day3).bi,99.99);
%     bi_percentile_high=mean([bi_percentile_high1 bi_percentile_high2  bi_percentile_high3]);
% 
%     for session = 1:24
%         % normalize muscle traces
%         tri_norm=(muscimol_cc_data.(mouse_num)(session).triM-tri_percentile_low)/(tri_percentile_high-tri_percentile_low);
%         tri_norm(tri_norm<0)=0;
%         muscimol_cc_data.(mouse_num)(session).tri_norm=tri_norm;
%         bi_norm=(muscimol_cc_data.(mouse_num)(session).bi-bi_percentile_low)/(bi_percentile_high-bi_percentile_low);
%         bi_norm(bi_norm<0)=0;
%         muscimol_cc_data.(mouse_num)(session).bi_norm=bi_norm;
%         % calculate CCI
%         CCI=nan*ones(1,length(muscimol_cc_data.(mouse_num)(session).bi_norm));
%         for k=1:length(CCI)
%             if muscimol_cc_data.(mouse_num)(session).bi_norm(k) > muscimol_cc_data.(mouse_num)(session).tri_norm(k)
%                 CCI(k)= (muscimol_cc_data.(mouse_num)(session).tri_norm(k)/muscimol_cc_data.(mouse_num)(session).bi_norm(k))*...
%                 (muscimol_cc_data.(mouse_num)(session).bi_norm(k)+muscimol_cc_data.(mouse_num)(session).tri_norm(k));
%             elseif muscimol_cc_data.(mouse_num)(session).tri_norm(k) > muscimol_cc_data.(mouse_num)(session).bi_norm(k)    
%             CCI(k)= (muscimol_cc_data.(mouse_num)(session).bi_norm(k)/muscimol_cc_data.(mouse_num)(session).tri_norm(k))*...
%                 (muscimol_cc_data.(mouse_num)(session).bi_norm(k)+muscimol_cc_data.(mouse_num)(session).tri_norm(k));
%             else %if they're the same, ratio will be 1
%                 CCI(k)= (muscimol_cc_data.(mouse_num)(session).bi_norm(k)/muscimol_cc_data.(mouse_num)(session).tri_norm(k))*...
%                 (muscimol_cc_data.(mouse_num)(session).bi_norm(k)+muscimol_cc_data.(mouse_num)(session).tri_norm(k));
%             end
%         end
%         muscimol_cc_data.(mouse_num)(session).CCI=CCI;
%         muscimol_cc_data.(mouse_num)(session).CCI_mean=mean(CCI);
%         muscimol_cc_data.(mouse_num)(session).CCI_median=median(CCI);
%         muscimol_cc_data.(mouse_num)(session).CCI_pmax=prctile(CCI,99.99); 
%     end
% end
% 
% % since only A37 has a 24th session, have it replace the 23 and delete the empty row 24 from others
% muscimol_cc_data.A37(23)=muscimol_cc_data.A37(24);
% muscimol_cc_data.A37(24)=[];
% muscimol_cc_data.A38(24)=[];
% muscimol_cc_data.A39(24)=[];
% muscimol_cc_data.A40(24)=[];
% muscimol_cc_data.A41(24)=[];
% 
% %save data
% save('muscimol_cc_data','muscimol_cc_data','-v7.3');
% cd('F:');
% 

%% start here if you opened prev saved file


%%  COCON LEARNING (muscimol cohort)

% alt learning metric already in the matrix

% cc learning, notes on parameter choice:
% the perc_thr_onset of 0.5 matches the alignment method of the trial selection (but that uses triceps trace instead of CCI)
% the qui_inds match the trial selection ones: 100 ms before onset, quiescent for at least 50 ms
% peak theshold for CCI event is 0.2 (no need for equivalency here since different metric than pure EMG)
% CCI quiencence threshold is 0.1

window=500;
perc_thr_onset=0.5;
peak_th=0.2;
qui_thresh_cci=0.1;
qui_inds=350:400; 
qui_thresh=0.25;
mouse_IDs={'A37','A38','A39','A40','A41'};

for i=1:length(mouse_IDs)
    i
    mouse_num=mouse_IDs{i};

     for session = 1:23
         session
        % look at number of CCI events datapoints above a threshold
        [pks,locs]=findpeaks(muscimol_cc_data.(mouse_num)(session).CCI,'MinPeakHeight',0.1,'MinPeakDistance',800,'MinPeakProminence',peak_th);
        pks((locs-window)<0)=[];
        locs((locs-window)<0)=[]; % remove those locs for whom locs-window<0
        CCI_peaks_mat=nan*ones(length(locs),window);
        for ii=1:length(locs)
            CCI_peaks_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).CCI(locs(ii)-window+1:locs(ii));
        end
        num_events_no_qui=length(locs);
        
        % align on rising edge
        locs_align=[];
        for ii=1:length(locs)
            trial=CCI_peaks_mat(ii,:);
            peak_height=pks(ii);
            idx=find(trial<perc_thr_onset*peak_height,1,'last');
            if ~isempty(idx)
                locs_align=[locs_align (locs(ii)-(window-idx))];
            end
        end
        locs_align((locs_align-window)<0)=[]; % remove those locs for whom locs-window<0

        %remove peaks that have non-quiescent preceding activity in EMG AND remove those with preceding cocon
        tri_mat=nan*ones(length(locs_align),window);
        bi_mat=nan*ones(length(locs_align),window);
        CCI_peaks_mat=nan*ones(length(locs_align),window);
        for ii=1:length(locs_align)
            tri_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).tri_norm(locs_align(ii)-window+1:locs_align(ii));
            bi_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).bi_norm(locs_align(ii)-window+1:locs_align(ii));
            CCI_peaks_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).CCI(locs_align(ii)-window+1:locs_align(ii));
        end
        
        throw_vec=[];
        for ii=1:length(locs_align)
            if mean(tri_mat(ii,qui_inds))>qui_thresh || mean(bi_mat(ii,qui_inds))>qui_thresh
                throw_vec=[throw_vec ii];
            end
            if mean(CCI_peaks_mat(ii,qui_inds))>qui_thresh_cci
                throw_vec=[throw_vec ii];
            end
        end
        throw_vec=unique(throw_vec);
        locs_align(throw_vec)=[];

        muscimol_cc_data.(mouse_num)(session).CCI_events=length(locs_align);
        muscimol_cc_data.(mouse_num)(session).CCI_event_indx=locs_align;

%         % plot the events
%         tri_mat=nan*ones(length(locs_align),window*2);
%         bi_mat=nan*ones(length(locs_align),window*2);
%         %CCI_peaks_mat=nan*ones(length(locs_align),window);
%         for ii=1:length(locs_align)
%             tri_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).tri_norm(locs_align(ii)-window+1:locs_align(ii)+window);
%             bi_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).bi_norm(locs_align(ii)-window+1:locs_align(ii)+window);
%             %CCI_peaks_mat(ii,:)=muscimol_cc_data.(mouse_num)(session).CCI(locs_align(ii)-window+1:locs_align(ii));
%         end
%         
%         subplot(4,6,session);
%         hold on;
%         shadedErrorBar(1:window*2,mean(tri_mat),std(tri_mat)/sqrt(size(tri_mat,1)),{'Color',[0 1 0],'LineWidth',2},1);
%         shadedErrorBar(1:window*2,mean(bi_mat),std(bi_mat)/sqrt(size(bi_mat,1)),{'Color',[0 0 1],'LineWidth',2},1);
%         axis square;
%         xlim([250 750]);
    end
end


%% learning curves

% alt
param='velMean20min';
mus_cohort_mat=[vertcat(muscimol_alt_data.A37(1:14).(param))'; vertcat(muscimol_alt_data.A38(1:14).(param))'; vertcat(muscimol_alt_data.A39(1:14).(param))';... 
    vertcat(muscimol_alt_data.A40(1:14).(param))'; vertcat(muscimol_alt_data.A41(1:14).(param))'];
mus_cohort_mat=[mus_cohort_mat nan*ones(5,1);];

% cc
learning_mat=[vertcat(muscimol_cc_data.A37.CCI_events) vertcat(muscimol_cc_data.A38.CCI_events) ...
    vertcat(muscimol_cc_data.A39.CCI_events) vertcat(muscimol_cc_data.A40.CCI_events) ...
    vertcat(muscimol_cc_data.A41.CCI_events)];

figure;
%alt traces
subplot(2,4,1); hold on;
plot(mus_cohort_mat','Color',[0.5 0.5 0.5]);
plot(mean(mus_cohort_mat),'b','LineWidth',2);
ylabel('Mean wheel velocity (cm/s)');
ylim([0 3.1]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;
title('Single traces');

% sem alt mean
subplot(2,4,2); hold on;
shadedErrorBar(1:15,mean(mus_cohort_mat),std(mus_cohort_mat)/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Mean wheel velocity (cm/s)');
ylim([0 2.5]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;
title('Mean with SEM');

%nomed to D1 performance
mus_cohort_mat_norm=mus_cohort_mat./mus_cohort_mat(:,1);
subplot(2,4,3); hold on;
shadedErrorBar(1:15,mean(mus_cohort_mat_norm),std(mus_cohort_mat_norm)/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Norm mean velocity (cm/s)');
ylim([1 20]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;
title('Normed to D1');

%nomed to mean D1-2 performance
mus_cohort_mat_norm=mus_cohort_mat./mean([mus_cohort_mat(:,1) mus_cohort_mat(:,2)],2);
subplot(2,4,4); hold on;
shadedErrorBar(1:15,mean(mus_cohort_mat_norm),std(mus_cohort_mat_norm)/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Norm mean velocity (cm/s)');
ylim([0.4 8]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;
title('Normed to D1-D2 mean');

% cocon
subplot(2,4,5); hold on;
plot(learning_mat,'Color',[0.5 0.5 0.5]);
plot(mean(learning_mat,2),'b','LineWidth',2);
ylabel('Number of cocon events');
ylim([0 252]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;

    
subplot(2,4,6); hold on;
shadedErrorBar(1:23,mean(learning_mat'),std(learning_mat')/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Number of cocon events');
ylim([0 160]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;

% normed to D1 performance
learning_mat_norm=learning_mat./learning_mat(1,:);

subplot(2,4,7); hold on;
shadedErrorBar(1:23,mean(learning_mat_norm'),std(learning_mat_norm')/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Norm num cocon events');
ylim([0.7 2.02]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;

% normed to mean D1-2 performance
learning_mat_norm=learning_mat./mean([learning_mat(1,:); learning_mat(2,:)]);
subplot(2,4,8); hold on;
shadedErrorBar(1:23,mean(learning_mat_norm'),std(learning_mat_norm')/sqrt(5),{'Color',[0 0 0],'LineWidth',2},1);
ylabel('Norm num cocon events');
ylim([0.8 2.3]);

xlim([0 15]);
xticks([1 5 10 15]);
xlabel('Session number');
axis square;
box off;

suptitle('learning curves from muscimol cohort');

%% plot muscimol quantification

% alt
param='velMean20min';
baseline_inds=12:14; % 3 days prior
show_inds=15:22;
timecourse=[vertcat(muscimol_alt_data.A37(show_inds).(param))/mean(vertcat(muscimol_alt_data.A37(baseline_inds).(param))) ...
    vertcat(muscimol_alt_data.A38(show_inds).(param))/mean(vertcat(muscimol_alt_data.A38(baseline_inds).(param))) ...
    vertcat(muscimol_alt_data.A39(show_inds).(param))/mean(vertcat(muscimol_alt_data.A39(baseline_inds).(param))) ...
    vertcat(muscimol_alt_data.A40(show_inds).(param))/mean(vertcat(muscimol_alt_data.A40(baseline_inds).(param))) ...
    vertcat(muscimol_alt_data.A41(show_inds).(param))/mean(vertcat(muscimol_alt_data.A41(baseline_inds).(param)))];
contra_musc_alt=[timecourse(1,1); timecourse(1,2); timecourse(1,3); timecourse(3,4); timecourse(3,5)];
contra_musc_wash_alt=[timecourse(2,1); timecourse(2,2); timecourse(2,3); timecourse(4,4); timecourse(4,5)];
contra_sal_alt=[timecourse(3,1); timecourse(3,2); timecourse(3,3); timecourse(1,4); timecourse(1,5)];
contra_sal_wash_alt=[timecourse(4,1); timecourse(4,2); timecourse(4,3); timecourse(2,4); timecourse(2,5)];
bilat_musc_alt=[timecourse(5,1); timecourse(5,2); timecourse(5,3); timecourse(7,4); timecourse(7,5)];
bilat_musc_wash_alt=[timecourse(6,1); timecourse(6,2); timecourse(6,3); timecourse(8,4); timecourse(8,5)];
bilat_sal_alt=[timecourse(7,1); timecourse(7,2); timecourse(7,3); timecourse(5,4); timecourse(5,5)];
bilat_sal_wash_alt=[timecourse(8,1); timecourse(8,2); timecourse(8,3); timecourse(6,4); timecourse(6,5)];

figure;
%just contralateral
subplot(2,3,1); 
hold on;
errorbar(1,mean(contra_musc_alt),std(contra_musc_alt)/sqrt(5),'k');
plot(1,mean(contra_musc_alt),'.k');
errorbar(2,mean(contra_musc_wash_alt),std(contra_musc_wash_alt)/sqrt(5),'k');
plot(2,mean(contra_musc_wash_alt),'.k');
errorbar(3,mean(contra_sal_alt),std(contra_sal_alt)/sqrt(5),'k');
plot(3,mean(contra_sal_alt),'.k');
errorbar(4,mean(contra_sal_wash_alt),std(contra_sal_wash_alt)/sqrt(5),'k');
plot(4,mean(contra_sal_wash_alt),'.k');
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_alt) mean(contra_musc_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([3 4],[mean(contra_sal_alt) mean(contra_sal_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([0.5 4.5],[1 1],'--k');
xlim([0.5 4.5]);
ylim([0 1.5]);
axis square;
box off;
xticks([1:8]);
xticklabels({'mus contra','washout','sal contra','washout','mus bilat','washout','sal bilat','washout'});
xtickangle(45); 
ylabel('Wheel velocity/baseline');
title('Alt, contra');

%contralateral with bilat overlay
gray=[0.7 0.7 0.7];
subplot(2,3,2); 
hold on;
errorbar(1,mean(contra_musc_alt),std(contra_musc_alt)/sqrt(5),'Color',gray);
plot(1,mean(contra_musc_alt),'.','Color',gray);
errorbar(2,mean(contra_musc_wash_alt),std(contra_musc_wash_alt)/sqrt(5),'Color',gray);
plot(2,mean(contra_musc_wash_alt),'.','Color',gray);
errorbar(3,mean(contra_sal_alt),std(contra_sal_alt)/sqrt(5),'Color',gray);
plot(3,mean(contra_sal_alt),'.','Color',gray);
errorbar(4,mean(contra_sal_wash_alt),std(contra_sal_wash_alt)/sqrt(5),'Color',gray);
plot(4,mean(contra_sal_wash_alt),'.','Color',gray);
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_alt) mean(contra_musc_wash_alt)],'Color',gray);
plot([3 4],[mean(contra_sal_alt) mean(contra_sal_wash_alt)],'Color',gray);
% bilat overlay
errorbar(1,mean(bilat_musc_alt),std(bilat_musc_alt)/sqrt(5),'k');
plot(1,mean(bilat_musc_alt),'.k');
errorbar(2,mean(bilat_musc_wash_alt),std(bilat_musc_wash_alt)/sqrt(5),'k');
plot(2,mean(bilat_musc_wash_alt),'.k');
errorbar(3,mean(bilat_sal_alt),std(bilat_sal_alt)/sqrt(5),'k');
plot(3,mean(bilat_sal_alt),'.k');
errorbar(4,mean(bilat_sal_wash_alt),std(bilat_sal_wash_alt)/sqrt(5),'k');
plot(4,mean(bilat_sal_wash_alt),'.k');
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(bilat_musc_alt) mean(bilat_musc_wash_alt)],'k');
plot([3 4],[mean(bilat_sal_alt) mean(bilat_sal_wash_alt)],'k');
plot([0.5 4.5],[1 1],'--k');
xlim([0.5 4.5]);
ylim([0 1.5]);
axis square;
box off;
xticks([1:8]);
xticklabels({'mus contra','washout','sal contra','washout','mus bilat','washout','sal bilat','washout'});
xtickangle(45); 
ylabel('Wheel velocity/baseline');
title('Contra + bilat overlay');


% contra and bilateral side by side
subplot(2,3,3); 
hold on;
errorbar(1,mean(contra_musc_alt),std(contra_musc_alt)/sqrt(5),'k');
plot(1,mean(contra_musc_alt),'.k');
errorbar(2,mean(contra_musc_wash_alt),std(contra_musc_wash_alt)/sqrt(5),'k');
plot(2,mean(contra_musc_wash_alt),'.k');
errorbar(3,mean(contra_sal_alt),std(contra_sal_alt)/sqrt(5),'k');
plot(3,mean(contra_sal_alt),'.k');
errorbar(4,mean(contra_sal_wash_alt),std(contra_sal_wash_alt)/sqrt(5),'k');
plot(4,mean(contra_sal_wash_alt),'.k');
errorbar(5,mean(bilat_musc_alt),std(bilat_musc_alt)/sqrt(5),'k');
plot(5,mean(bilat_musc_alt),'.k');
errorbar(6,mean(bilat_musc_wash_alt),std(bilat_musc_wash_alt)/sqrt(5),'k');
plot(6,mean(bilat_musc_wash_alt),'.k');
errorbar(7,mean(bilat_sal_alt),std(bilat_sal_alt)/sqrt(5),'k');
plot(7,mean(bilat_sal_alt),'.k');
errorbar(8,mean(bilat_sal_wash_alt),std(bilat_sal_wash_alt)/sqrt(5),'k');
plot(8,mean(bilat_sal_wash_alt),'.k');
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_alt) mean(contra_musc_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([3 4],[mean(contra_sal_alt) mean(contra_sal_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([5 6],[mean(bilat_musc_alt) mean(bilat_musc_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([7 8],[mean(bilat_sal_alt) mean(bilat_sal_wash_alt)],'Color',[0.5 0.5 0.5]);
plot([0.5 8.5],[1 1],'--k');
xlim([0.5 8.5]);
ylim([0 1.5]);
axis square;
box off;
xticks([1:8]);
xticklabels({'mus contra','washout','sal contra','washout','mus bilat','washout','sal bilat','washout'});
xtickangle(45); 
ylabel('Wheel velocity/baseline');
title('Contra and bilat');



[h,p]=ttest(contra_musc_alt,contra_sal_alt,'Tail','left'); %0.0228
[h,p]=ttest(bilat_musc_alt,bilat_sal_alt,'Tail','left'); %0.0024

[h,p]=ttest(contra_musc_alt,contra_musc_wash_alt,'Tail','left'); %0.0017
[h,p]=ttest(contra_sal_alt,contra_sal_wash_alt,'Tail','left'); %0.6424

[h,p]=ttest(bilat_musc_alt,bilat_musc_wash_alt,'Tail','left'); %0.0034
[h,p]=ttest(bilat_sal_alt,bilat_sal_wash_alt,'Tail','left'); %0.2249

% cocon
baseline_inds=13:15; % 3 days prior
show_inds=16:23;
gray=[0.7 0.7 0.7];

timecourse=[vertcat(muscimol_cc_data.A37([16:23]).CCI_events)/mean(vertcat(muscimol_cc_data.A37(baseline_inds).CCI_events)) ...
    vertcat(muscimol_cc_data.A38(show_inds).CCI_events)/mean(vertcat(muscimol_cc_data.A38(baseline_inds).CCI_events)) ...
    vertcat(muscimol_cc_data.A39(show_inds).CCI_events)/mean(vertcat(muscimol_cc_data.A39(baseline_inds).CCI_events)) ...
    vertcat(muscimol_cc_data.A40(show_inds).CCI_events)/mean(vertcat(muscimol_cc_data.A40(baseline_inds).CCI_events)) ...
    vertcat(muscimol_cc_data.A41(show_inds).CCI_events)/mean(vertcat(muscimol_cc_data.A41(baseline_inds).CCI_events))];
contra_musc_cc=[timecourse(1,1); timecourse(1,2); timecourse(1,3); timecourse(3,4); timecourse(3,5)];
contra_musc_wash_cc=[timecourse(2,1); timecourse(2,2); timecourse(2,3); timecourse(4,4); timecourse(4,5)];
contra_sal_cc=[timecourse(3,1); timecourse(3,2); timecourse(3,3); timecourse(1,4); timecourse(1,5)];
contra_sal_wash_cc=[timecourse(4,1); timecourse(4,2); timecourse(4,3); timecourse(2,4); timecourse(2,5)];
bilat_musc_cc=[timecourse(5,1); timecourse(5,2); timecourse(5,3); timecourse(7,4); timecourse(7,5)];
bilat_musc_wash_cc=[timecourse(6,1); timecourse(6,2); timecourse(6,3); timecourse(8,4); timecourse(8,5)];
bilat_sal_cc=[timecourse(7,1); timecourse(7,2); timecourse(7,3); timecourse(5,4); timecourse(5,5)];
bilat_sal_wash_cc=[timecourse(8,1); timecourse(8,2); timecourse(8,3); timecourse(6,4); timecourse(6,5)];

% contralateral
figure; hold on;
subplot(2,3,4); hold on;
errorbar(1,mean(contra_musc_cc),std(contra_musc_cc)/sqrt(5),'k');
plot(1,mean(contra_musc_cc),'.k');
errorbar(2,mean(contra_musc_wash_cc),std(contra_musc_wash_cc)/sqrt(5),'k');
plot(2,mean(contra_musc_wash_cc),'.k');
errorbar(3,mean(contra_sal_cc),std(contra_sal_cc)/sqrt(5),'k');
plot(3,mean(contra_sal_cc),'.k');
errorbar(4,mean(contra_sal_wash_cc),std(contra_sal_wash_cc)/sqrt(5),'k');
plot(4,mean(contra_sal_wash_cc),'.k');
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_cc) mean(contra_musc_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([3 4],[mean(contra_sal_cc) mean(contra_sal_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([0.5 4.5],[1 1],'--k');
ylim([0 1.5]);
xlim([0.5 4.5]);
axis square;
box off;
xticks([1:4]);
xticklabels({'mus contra','washout','sal contra','washout'});
xtickangle(45); 
ylabel('Number of cocons/baseline');
title('Cocon, contra');

% contralateral with bilat overlay
subplot(2,3,5); hold on;

errorbar(1,mean(contra_musc_cc),std(contra_musc_cc)/sqrt(5),'Color',gray);
plot(1,mean(contra_musc_cc),'.','Color',gray);
errorbar(2,mean(contra_musc_wash_cc),std(contra_musc_wash_cc)/sqrt(5),'Color',gray);
plot(2,mean(contra_musc_wash_cc),'.','Color',gray);
errorbar(3,mean(contra_sal_cc),std(contra_sal_cc)/sqrt(5),'Color',gray);
plot(3,mean(contra_sal_cc),'.','Color',gray);
errorbar(4,mean(contra_sal_wash_cc),std(contra_sal_wash_cc)/sqrt(5),'Color',gray);
plot(4,mean(contra_sal_wash_cc),'.','Color',gray);
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_cc) mean(contra_musc_wash_cc)],'Color',gray);
plot([3 4],[mean(contra_sal_cc) mean(contra_sal_wash_cc)],'Color',gray);
% overlay bilat
errorbar(1,mean(bilat_musc_cc),std(bilat_musc_cc)/sqrt(5),'k');
plot(1,mean(bilat_musc_cc),'.k');
errorbar(2,mean(bilat_musc_wash_cc),std(bilat_musc_wash_cc)/sqrt(5),'k');
plot(2,mean(bilat_musc_wash_cc),'.k');
errorbar(3,mean(bilat_sal_cc),std(bilat_sal_cc)/sqrt(5),'k');
plot(3,mean(bilat_sal_cc),'.k');
errorbar(4,mean(bilat_sal_wash_cc),std(bilat_sal_wash_cc)/sqrt(5),'k');
plot(4,mean(bilat_sal_wash_cc),'.k');
plot([1 2],[mean(bilat_musc_cc) mean(bilat_musc_wash_cc)],'k');
plot([3 4],[mean(bilat_sal_cc) mean(bilat_sal_wash_cc)],'k');

plot([0.5 4.5],[1 1],'--k');
ylim([0 2]);
xlim([0.5 4.5]);
axis square;
box off;
xticks([1:4]);
xticklabels({'mus contra','washout','sal contra','washout'});
xtickangle(45); 
ylabel('Number of cocons/baseline');
title('Contra + bilat overlay');


subplot(2,3,6); hold on;
errorbar(1,mean(contra_musc_cc),std(contra_musc_cc)/sqrt(5),'k');
plot(1,mean(contra_musc_cc),'.k');
errorbar(2,mean(contra_musc_wash_cc),std(contra_musc_wash_cc)/sqrt(5),'k');
plot(2,mean(contra_musc_wash_cc),'.k');
errorbar(3,mean(contra_sal_cc),std(contra_sal_cc)/sqrt(5),'k');
plot(3,mean(contra_sal_cc),'.k');
errorbar(4,mean(contra_sal_wash_cc),std(contra_sal_wash_cc)/sqrt(5),'k');
plot(4,mean(contra_sal_wash_cc),'.k');
errorbar(5,mean(bilat_musc_cc),std(bilat_musc_cc)/sqrt(5),'k');
plot(5,mean(bilat_musc_cc),'.k');
errorbar(6,mean(bilat_musc_wash_cc),std(bilat_musc_wash_cc)/sqrt(5),'k');
plot(6,mean(bilat_musc_wash_cc),'.k');
errorbar(7,mean(bilat_sal_cc),std(bilat_sal_cc)/sqrt(5),'k');
plot(7,mean(bilat_sal_cc),'.k');
errorbar(8,mean(bilat_sal_wash_cc),std(bilat_sal_wash_cc)/sqrt(5),'k');
plot(8,mean(bilat_sal_wash_cc),'.k');
% with lines connecting musc-musc washout and sal-sal washout
plot([1 2],[mean(contra_musc_cc) mean(contra_musc_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([3 4],[mean(contra_sal_cc) mean(contra_sal_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([5 6],[mean(bilat_musc_cc) mean(bilat_musc_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([7 8],[mean(bilat_sal_cc) mean(bilat_sal_wash_cc)],'Color',[0.5 0.5 0.5]);
plot([0.5 8.5],[1 1],'--k');
xlim([0.5 8.5]);
ylim([0 2]);
axis square;
box off;
xticks([1:8]);
xticklabels({'mus contra','washout','sal contra','washout','mus bilat','washout','sal bilat','washout'});
xtickangle(45); 
ylabel('Number of cocons/baseline');
title('Contra and bilat');


[h,p]=ttest(contra_musc_cc,contra_sal_cc,'Tail','left'); %0.0042
[h,p]=ttest(bilat_musc_cc,bilat_sal_cc,'Tail','left'); %0.0673

[h,p]=ttest(contra_musc_cc,contra_musc_wash_cc,'Tail','left'); %0.0037
[h,p]=ttest(contra_sal_cc,contra_sal_wash_cc,'Tail','left'); %0.7656

[h,p]=ttest(bilat_musc_cc,bilat_musc_wash_cc,'Tail','left'); %0.0280
[h,p]=ttest(bilat_sal_cc,bilat_sal_wash_cc,'Tail','left'); %0.3499

% without outlier, so 4/5 showed significant deficit
[h,p]=ttest(bilat_musc_cc(1:4),bilat_sal_cc(1:4),'Tail','left'); %0.0154

suptitle('Effects of muscimol on performance');

%% deep dive on Alt

% when wheel moving FORWARD

%figure out threshold for wheel velocity, find baseline (no movement, and
%take the STD, then 5-7 STDs above that is your thresh, be conservative
inds_a=317000:338000;
inds_b=844000:857000;
inds_c=1090500:1098000;
inds_d=56000:67000;

% figure;
% plot(muscimol_alt_data.A37(13).enc_vel,'k');
% hold on;
% plot(inds_a,0.2*ones(1,length(inds_a)),'g');
% plot(inds_b,0.2*ones(1,length(inds_b)),'g');
% plot(inds_c,0.2*ones(1,length(inds_c)),'g');

low_trace=[muscimol_alt_data.A37(13).enc_vel(inds_a) muscimol_alt_data.A37(13).enc_vel(inds_b) ...
    muscimol_alt_data.A37(13).enc_vel(inds_c) muscimol_alt_data.A37(13).enc_vel(inds_d)];
std_low=std(low_trace);
th=std_low*5; %7


% figure;
% plot(muscimol_alt_data.A37(13).enc_vel,'k');
% hold on;
% plot(1:length(muscimol_alt_data.A37(13).enc_vel),th*ones(1,length(muscimol_alt_data.A37(13).enc_vel)),'g');

data_vel_bases=struct;
data_vel_cond=struct;
sesh_dur=1200000;

%alt baseline is sessions 12, 13, 14
%baseline session 1 (# 12)
base_sesh_num=12;
data_vel_bases(1).baseline1=muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);
data_vel_bases(2).baseline1=muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);
data_vel_bases(3).baseline1=muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);
data_vel_bases(4).baseline1=muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);
data_vel_bases(5).baseline1=muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);

data_vel_bases(1).baseline1_mean=mean(muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).baseline1_mean=mean(muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).baseline1_mean=mean(muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).baseline1_mean=mean(muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).baseline1_mean=mean(muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).baseline1_perc=length(data_vel_bases(1).baseline1)/sesh_dur;
data_vel_bases(2).baseline1_perc=length(data_vel_bases(2).baseline1)/sesh_dur;
data_vel_bases(3).baseline1_perc=length(data_vel_bases(3).baseline1)/sesh_dur;
data_vel_bases(4).baseline1_perc=length(data_vel_bases(4).baseline1)/sesh_dur;
data_vel_bases(5).baseline1_perc=length(data_vel_bases(5).baseline1)/sesh_dur;

%baseline session 2 (# 13)
base_sesh_num=13;
data_vel_bases(1).baseline2=muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);
data_vel_bases(2).baseline2=muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);
data_vel_bases(3).baseline2=muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);
data_vel_bases(4).baseline2=muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);
data_vel_bases(5).baseline2=muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);

data_vel_bases(1).baseline2_mean=mean(muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).baseline2_mean=mean(muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).baseline2_mean=mean(muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).baseline2_mean=mean(muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).baseline2_mean=mean(muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).baseline2_perc=length(data_vel_bases(1).baseline2)/sesh_dur;
data_vel_bases(2).baseline2_perc=length(data_vel_bases(2).baseline2)/sesh_dur;
data_vel_bases(3).baseline2_perc=length(data_vel_bases(3).baseline2)/sesh_dur;
data_vel_bases(4).baseline2_perc=length(data_vel_bases(4).baseline2)/sesh_dur;
data_vel_bases(5).baseline2_perc=length(data_vel_bases(5).baseline2)/sesh_dur;

%baseline session 3 (# 14)
base_sesh_num=14;
data_vel_bases(1).baseline3=muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);
data_vel_bases(2).baseline3=muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);
data_vel_bases(3).baseline3=muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);
data_vel_bases(4).baseline3=muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);
data_vel_bases(5).baseline3=muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);

data_vel_bases(1).baseline3_mean=mean(muscimol_alt_data.A37(base_sesh_num).enc_vel(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).baseline3_mean=mean(muscimol_alt_data.A38(base_sesh_num).enc_vel(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).baseline3_mean=mean(muscimol_alt_data.A39(base_sesh_num).enc_vel(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).baseline3_mean=mean(muscimol_alt_data.A40(base_sesh_num).enc_vel(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).baseline3_mean=mean(muscimol_alt_data.A41(base_sesh_num).enc_vel(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).baseline3_perc=length(data_vel_bases(1).baseline3)/sesh_dur;
data_vel_bases(2).baseline3_perc=length(data_vel_bases(2).baseline3)/sesh_dur;
data_vel_bases(3).baseline3_perc=length(data_vel_bases(3).baseline3)/sesh_dur;
data_vel_bases(4).baseline3_perc=length(data_vel_bases(4).baseline3)/sesh_dur;
data_vel_bases(5).baseline3_perc=length(data_vel_bases(5).baseline3)/sesh_dur;


% average across the 3 baseline sessions
% data_vel_cond(1).baseline=mean([data_vel_bases(1)])
% data_vel_cond(2).baseline=muscimol_alt_data.A38(14).enc_vel(muscimol_alt_data.A38(14).enc_vel > th);
% data_vel_cond(3).baseline=muscimol_alt_data.A39(14).enc_vel(muscimol_alt_data.A39(14).enc_vel > th);
% data_vel_cond(4).baseline=muscimol_alt_data.A40(14).enc_vel(muscimol_alt_data.A40(14).enc_vel > th);
% data_vel_cond(5).baseline=muscimol_alt_data.A41(14).enc_vel(muscimol_alt_data.A41(14).enc_vel > th);

data_vel_cond(1).baseline_mean=mean([data_vel_bases(1).baseline1_mean data_vel_bases(1).baseline2_mean data_vel_bases(1).baseline3_mean]);
data_vel_cond(2).baseline_mean=mean([data_vel_bases(2).baseline1_mean data_vel_bases(2).baseline2_mean data_vel_bases(2).baseline3_mean]);
data_vel_cond(3).baseline_mean=mean([data_vel_bases(3).baseline1_mean data_vel_bases(3).baseline2_mean data_vel_bases(3).baseline3_mean]);
data_vel_cond(4).baseline_mean=mean([data_vel_bases(4).baseline1_mean data_vel_bases(4).baseline2_mean data_vel_bases(4).baseline3_mean]);
data_vel_cond(5).baseline_mean=mean([data_vel_bases(5).baseline1_mean data_vel_bases(5).baseline2_mean data_vel_bases(5).baseline3_mean]);

data_vel_cond(1).baseline_perc=mean([data_vel_bases(1).baseline1_perc data_vel_bases(1).baseline2_perc data_vel_bases(1).baseline3_perc]);
data_vel_cond(2).baseline_perc=mean([data_vel_bases(2).baseline1_perc data_vel_bases(2).baseline2_perc data_vel_bases(2).baseline3_perc]);
data_vel_cond(3).baseline_perc=mean([data_vel_bases(3).baseline1_perc data_vel_bases(3).baseline2_perc data_vel_bases(3).baseline3_perc]);
data_vel_cond(4).baseline_perc=mean([data_vel_bases(4).baseline1_perc data_vel_bases(4).baseline2_perc data_vel_bases(4).baseline3_perc]);
data_vel_cond(5).baseline_perc=mean([data_vel_bases(5).baseline1_perc data_vel_bases(5).baseline2_perc data_vel_bases(5).baseline3_perc]);



data_vel_cond(1).saline=muscimol_alt_data.A37(17).enc_vel(muscimol_alt_data.A37(17).enc_vel > th);
data_vel_cond(2).saline=muscimol_alt_data.A38(17).enc_vel(muscimol_alt_data.A38(17).enc_vel > th);
data_vel_cond(3).saline=muscimol_alt_data.A39(17).enc_vel(muscimol_alt_data.A39(17).enc_vel > th);
data_vel_cond(4).saline=muscimol_alt_data.A40(15).enc_vel(muscimol_alt_data.A40(15).enc_vel > th);
data_vel_cond(5).saline=muscimol_alt_data.A41(15).enc_vel(muscimol_alt_data.A41(15).enc_vel > th);

data_vel_cond(1).saline_mean=mean(muscimol_alt_data.A37(17).enc_vel(muscimol_alt_data.A37(17).enc_vel > th));
data_vel_cond(2).saline_mean=mean(muscimol_alt_data.A38(17).enc_vel(muscimol_alt_data.A38(17).enc_vel > th));
data_vel_cond(3).saline_mean=mean(muscimol_alt_data.A39(17).enc_vel(muscimol_alt_data.A39(17).enc_vel > th));
data_vel_cond(4).saline_mean=mean(muscimol_alt_data.A40(15).enc_vel(muscimol_alt_data.A40(15).enc_vel > th));
data_vel_cond(5).saline_mean=mean(muscimol_alt_data.A41(15).enc_vel(muscimol_alt_data.A41(15).enc_vel > th));

data_vel_cond(1).saline_perc=length(data_vel_cond(1).saline)/sesh_dur;
data_vel_cond(2).saline_perc=length(data_vel_cond(2).saline)/sesh_dur;
data_vel_cond(3).saline_perc=length(data_vel_cond(3).saline)/sesh_dur;
data_vel_cond(4).saline_perc=length(data_vel_cond(4).saline)/sesh_dur;
data_vel_cond(5).saline_perc=length(data_vel_cond(5).saline)/sesh_dur;

data_vel_cond(1).mus=muscimol_alt_data.A37(15).enc_vel(muscimol_alt_data.A37(15).enc_vel > th);
data_vel_cond(2).mus=muscimol_alt_data.A38(15).enc_vel(muscimol_alt_data.A38(15).enc_vel > th);
data_vel_cond(3).mus=muscimol_alt_data.A39(15).enc_vel(muscimol_alt_data.A39(15).enc_vel > th);
data_vel_cond(4).mus=muscimol_alt_data.A40(17).enc_vel(muscimol_alt_data.A40(17).enc_vel > th);
data_vel_cond(5).mus=muscimol_alt_data.A41(17).enc_vel(muscimol_alt_data.A41(17).enc_vel > th);

data_vel_cond(1).mus_mean=mean(muscimol_alt_data.A37(15).enc_vel(muscimol_alt_data.A37(15).enc_vel > th));
data_vel_cond(2).mus_mean=mean(muscimol_alt_data.A38(15).enc_vel(muscimol_alt_data.A38(15).enc_vel > th));
data_vel_cond(3).mus_mean=mean(muscimol_alt_data.A39(15).enc_vel(muscimol_alt_data.A39(15).enc_vel > th));
data_vel_cond(4).mus_mean=mean(muscimol_alt_data.A40(17).enc_vel(muscimol_alt_data.A40(17).enc_vel > th));
data_vel_cond(5).mus_mean=mean(muscimol_alt_data.A41(17).enc_vel(muscimol_alt_data.A41(17).enc_vel > th));

data_vel_cond(1).mus_perc=length(data_vel_cond(1).mus)/sesh_dur;
data_vel_cond(2).mus_perc=length(data_vel_cond(2).mus)/sesh_dur;
data_vel_cond(3).mus_perc=length(data_vel_cond(3).mus)/sesh_dur;
data_vel_cond(4).mus_perc=length(data_vel_cond(4).mus)/sesh_dur;
data_vel_cond(5).mus_perc=length(data_vel_cond(5).mus)/sesh_dur;


% tri and bi values when wheel engaged
% look across the 3 baseline sessions
% baseline session 1 (# 12)
% mean muscle activity
base_sesh_num=12;
data_vel_bases(1).base1_mean_tri=mean(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base1_mean_tri=mean(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base1_mean_tri=mean(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base1_mean_tri=mean(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base1_mean_tri=mean(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).base1_mean_bi=mean(muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base1_mean_bi=mean(muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base1_mean_bi=mean(muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base1_mean_bi=mean(muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base1_mean_bi=mean(muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri/bi pdist divided by length of data
data_vel_bases(1).base1_tribi_pdist=pdist([muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));

data_vel_bases(2).base1_tribi_pdist=pdist([muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));

data_vel_bases(3).base1_tribi_pdist=pdist([muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));

data_vel_bases(4).base1_tribi_pdist=pdist([muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
    
data_vel_bases(5).base1_tribi_pdist=pdist([muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri-bi correlation
corr_val=corrcoef(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(1).base1_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(2).base1_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(3).base1_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(4).base1_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));
data_vel_bases(5).base1_tribi_corr=corr_val(1,2);


%baseline session 2 (# 13)
base_sesh_num=13;
data_vel_bases(1).base2_mean_tri=mean(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base2_mean_tri=mean(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base2_mean_tri=mean(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base2_mean_tri=mean(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base2_mean_tri=mean(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).base2_mean_bi=mean(muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base2_mean_bi=mean(muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base2_mean_bi=mean(muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base2_mean_bi=mean(muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base2_mean_bi=mean(muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri/bi pdist
data_vel_bases(1).base2_tribi_pdist=pdist([muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));

data_vel_bases(2).base2_tribi_pdist=pdist([muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));

data_vel_bases(3).base2_tribi_pdist=pdist([muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));

data_vel_bases(4).base2_tribi_pdist=pdist([muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));

data_vel_bases(5).base2_tribi_pdist=pdist([muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri-bi correlation
corr_val=corrcoef(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(1).base2_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(2).base2_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(3).base2_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(4).base2_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));
data_vel_bases(5).base2_tribi_corr=corr_val(1,2);


%baseline session 3 (# 14)
base_sesh_num=14;
data_vel_bases(1).base3_mean_tri=mean(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base3_mean_tri=mean(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base3_mean_tri=mean(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base3_mean_tri=mean(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base3_mean_tri=mean(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

data_vel_bases(1).base3_mean_bi=mean(muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(2).base3_mean_bi=mean(muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(3).base3_mean_bi=mean(muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(4).base3_mean_bi=mean(muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(5).base3_mean_bi=mean(muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri/bi pdist
data_vel_bases(1).base3_tribi_pdist=pdist([muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));

data_vel_bases(2).base3_tribi_pdist=pdist([muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));

data_vel_bases(3).base3_tribi_pdist=pdist([muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));

data_vel_bases(4).base3_tribi_pdist=pdist([muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));

data_vel_bases(5).base3_tribi_pdist=pdist([muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th);...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));

% tri-bi correlation
corr_val=corrcoef(muscimol_alt_data.A37(base_sesh_num).tri_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A37(base_sesh_num).bi_norm(muscimol_alt_data.A37(base_sesh_num).enc_vel > th));
data_vel_bases(1).base3_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A38(base_sesh_num).tri_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A38(base_sesh_num).bi_norm(muscimol_alt_data.A38(base_sesh_num).enc_vel > th));
data_vel_bases(2).base3_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A39(base_sesh_num).tri_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A39(base_sesh_num).bi_norm(muscimol_alt_data.A39(base_sesh_num).enc_vel > th));
data_vel_bases(3).base3_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A40(base_sesh_num).tri_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A40(base_sesh_num).bi_norm(muscimol_alt_data.A40(base_sesh_num).enc_vel > th));
data_vel_bases(4).base3_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A41(base_sesh_num).tri_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th),...
    muscimol_alt_data.A41(base_sesh_num).bi_norm(muscimol_alt_data.A41(base_sesh_num).enc_vel > th));
data_vel_bases(5).base3_tribi_corr=corr_val(1,2);



% average over baselines
data_vel_cond(1).base_mean_tri=mean([data_vel_bases(1).base1_mean_tri data_vel_bases(1).base2_mean_tri data_vel_bases(1).base3_mean_tri]);
data_vel_cond(2).base_mean_tri=mean([data_vel_bases(2).base1_mean_tri data_vel_bases(2).base2_mean_tri data_vel_bases(2).base3_mean_tri]);
data_vel_cond(3).base_mean_tri=mean([data_vel_bases(3).base1_mean_tri data_vel_bases(3).base2_mean_tri data_vel_bases(3).base3_mean_tri]);
data_vel_cond(4).base_mean_tri=mean([data_vel_bases(4).base1_mean_tri data_vel_bases(4).base2_mean_tri data_vel_bases(4).base3_mean_tri]);
data_vel_cond(5).base_mean_tri=mean([data_vel_bases(5).base1_mean_tri data_vel_bases(5).base2_mean_tri data_vel_bases(5).base3_mean_tri]);

data_vel_cond(1).base_mean_bi=mean([data_vel_bases(1).base1_mean_bi data_vel_bases(1).base2_mean_bi data_vel_bases(1).base3_mean_bi]);
data_vel_cond(2).base_mean_bi=mean([data_vel_bases(2).base1_mean_bi data_vel_bases(2).base2_mean_bi data_vel_bases(2).base3_mean_bi]);
data_vel_cond(3).base_mean_bi=mean([data_vel_bases(3).base1_mean_bi data_vel_bases(3).base2_mean_bi data_vel_bases(3).base3_mean_bi]);
data_vel_cond(4).base_mean_bi=mean([data_vel_bases(4).base1_mean_bi data_vel_bases(4).base2_mean_bi data_vel_bases(4).base3_mean_bi]);
data_vel_cond(5).base_mean_bi=mean([data_vel_bases(5).base1_mean_bi data_vel_bases(5).base2_mean_bi data_vel_bases(5).base3_mean_bi]);

data_vel_cond(1).base_tribi_pdist=mean([data_vel_bases(1).base1_tribi_pdist data_vel_bases(1).base2_tribi_pdist data_vel_bases(1).base3_tribi_pdist]);
data_vel_cond(2).base_tribi_pdist=mean([data_vel_bases(2).base1_tribi_pdist data_vel_bases(2).base2_tribi_pdist data_vel_bases(2).base3_tribi_pdist]);
data_vel_cond(3).base_tribi_pdist=mean([data_vel_bases(3).base1_tribi_pdist data_vel_bases(3).base2_tribi_pdist data_vel_bases(3).base3_tribi_pdist]);
data_vel_cond(4).base_tribi_pdist=mean([data_vel_bases(4).base1_tribi_pdist data_vel_bases(4).base2_tribi_pdist data_vel_bases(4).base3_tribi_pdist]);
data_vel_cond(5).base_tribi_pdist=mean([data_vel_bases(5).base1_tribi_pdist data_vel_bases(5).base2_tribi_pdist data_vel_bases(5).base3_tribi_pdist]);

data_vel_cond(1).base_tribi_corr=mean([data_vel_bases(1).base1_tribi_corr data_vel_bases(1).base2_tribi_corr data_vel_bases(1).base3_tribi_corr]);
data_vel_cond(2).base_tribi_corr=mean([data_vel_bases(2).base1_tribi_corr data_vel_bases(2).base2_tribi_corr data_vel_bases(2).base3_tribi_corr]);
data_vel_cond(3).base_tribi_corr=mean([data_vel_bases(3).base1_tribi_corr data_vel_bases(3).base2_tribi_corr data_vel_bases(3).base3_tribi_corr]);
data_vel_cond(4).base_tribi_corr=mean([data_vel_bases(4).base1_tribi_corr data_vel_bases(4).base2_tribi_corr data_vel_bases(4).base3_tribi_corr]);
data_vel_cond(5).base_tribi_corr=mean([data_vel_bases(5).base1_tribi_corr data_vel_bases(5).base2_tribi_corr data_vel_bases(5).base3_tribi_corr]);



% sal
data_vel_cond(1).sal_mean_tri=mean(muscimol_alt_data.A37(17).tri_norm(muscimol_alt_data.A37(17).enc_vel > th));
data_vel_cond(2).sal_mean_tri=mean(muscimol_alt_data.A38(17).tri_norm(muscimol_alt_data.A38(17).enc_vel > th));
data_vel_cond(3).sal_mean_tri=mean(muscimol_alt_data.A39(17).tri_norm(muscimol_alt_data.A39(17).enc_vel > th));
data_vel_cond(4).sal_mean_tri=mean(muscimol_alt_data.A40(15).tri_norm(muscimol_alt_data.A40(15).enc_vel > th));
data_vel_cond(5).sal_mean_tri=mean(muscimol_alt_data.A41(15).tri_norm(muscimol_alt_data.A41(15).enc_vel > th));

data_vel_cond(1).sal_mean_bi=mean(muscimol_alt_data.A37(17).bi_norm(muscimol_alt_data.A37(17).enc_vel > th));
data_vel_cond(2).sal_mean_bi=mean(muscimol_alt_data.A38(17).bi_norm(muscimol_alt_data.A38(17).enc_vel > th));
data_vel_cond(3).sal_mean_bi=mean(muscimol_alt_data.A39(17).bi_norm(muscimol_alt_data.A39(17).enc_vel > th));
data_vel_cond(4).sal_mean_bi=mean(muscimol_alt_data.A40(15).bi_norm(muscimol_alt_data.A40(15).enc_vel > th));
data_vel_cond(5).sal_mean_bi=mean(muscimol_alt_data.A41(15).bi_norm(muscimol_alt_data.A41(15).enc_vel > th));

%pdist
sesh_num=17;
data_vel_cond(1).sal_tribi_pdist=pdist([muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th);...
    muscimol_alt_data.A37(sesh_num).bi_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th));
    
data_vel_cond(2).sal_tribi_pdist=pdist([muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th);...
    muscimol_alt_data.A38(sesh_num).bi_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th));
    
data_vel_cond(3).sal_tribi_pdist=pdist([muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th);...
    muscimol_alt_data.A39(sesh_num).bi_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th));
sesh_num=15;
data_vel_cond(4).sal_tribi_pdist=pdist([muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th);...
    muscimol_alt_data.A40(sesh_num).bi_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th));
    
data_vel_cond(5).sal_tribi_pdist=pdist([muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th);...
    muscimol_alt_data.A41(sesh_num).bi_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th));

% tri-bi correlation
sesh_num=17;
corr_val=corrcoef(muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th),...
    muscimol_alt_data.A37(sesh_num).bi_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th));
data_vel_cond(1).sal_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th),...
    muscimol_alt_data.A38(sesh_num).bi_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th));
data_vel_cond(2).sal_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th),...
    muscimol_alt_data.A39(sesh_num).bi_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th));
data_vel_cond(3).sal_tribi_corr=corr_val(1,2);

sesh_num=15;
corr_val=corrcoef(muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th),...
    muscimol_alt_data.A40(sesh_num).bi_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th));
data_vel_cond(4).sal_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th),...
    muscimol_alt_data.A41(sesh_num).bi_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th));
data_vel_cond(5).sal_tribi_corr=corr_val(1,2);


%mus
data_vel_cond(1).mus_mean_tri=mean(muscimol_alt_data.A37(15).tri_norm(muscimol_alt_data.A37(15).enc_vel > th));
data_vel_cond(2).mus_mean_tri=mean(muscimol_alt_data.A38(15).tri_norm(muscimol_alt_data.A38(15).enc_vel > th));
data_vel_cond(3).mus_mean_tri=mean(muscimol_alt_data.A39(15).tri_norm(muscimol_alt_data.A39(15).enc_vel > th));
data_vel_cond(4).mus_mean_tri=mean(muscimol_alt_data.A40(17).tri_norm(muscimol_alt_data.A40(17).enc_vel > th));
data_vel_cond(5).mus_mean_tri=mean(muscimol_alt_data.A41(17).tri_norm(muscimol_alt_data.A41(17).enc_vel > th));

data_vel_cond(1).mus_mean_bi=mean(muscimol_alt_data.A37(15).bi_norm(muscimol_alt_data.A37(15).enc_vel > th));
data_vel_cond(2).mus_mean_bi=mean(muscimol_alt_data.A38(15).bi_norm(muscimol_alt_data.A38(15).enc_vel > th));
data_vel_cond(3).mus_mean_bi=mean(muscimol_alt_data.A39(15).bi_norm(muscimol_alt_data.A39(15).enc_vel > th));
data_vel_cond(4).mus_mean_bi=mean(muscimol_alt_data.A40(17).bi_norm(muscimol_alt_data.A40(17).enc_vel > th));
data_vel_cond(5).mus_mean_bi=mean(muscimol_alt_data.A41(17).bi_norm(muscimol_alt_data.A41(17).enc_vel > th));

%pdist
sesh_num=15;
data_vel_cond(1).mus_tribi_pdist=pdist([muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th);...
    muscimol_alt_data.A37(sesh_num).bi_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th));
    
data_vel_cond(2).mus_tribi_pdist=pdist([muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th);...
    muscimol_alt_data.A38(sesh_num).bi_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th));
    
data_vel_cond(3).mus_tribi_pdist=pdist([muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th);...
    muscimol_alt_data.A39(sesh_num).bi_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th));
    
sesh_num=17;
data_vel_cond(4).mus_tribi_pdist=pdist([muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th);...
    muscimol_alt_data.A40(sesh_num).bi_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th));
    
data_vel_cond(5).mus_tribi_pdist=pdist([muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th);...
    muscimol_alt_data.A41(sesh_num).bi_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th)])...
    /length(muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th));

% tri-bi correlation
sesh_num=15;
corr_val=corrcoef(muscimol_alt_data.A37(sesh_num).tri_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th),...
    muscimol_alt_data.A37(sesh_num).bi_norm(muscimol_alt_data.A37(sesh_num).enc_vel > th));
data_vel_cond(1).mus_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A38(sesh_num).tri_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th),...
    muscimol_alt_data.A38(sesh_num).bi_norm(muscimol_alt_data.A38(sesh_num).enc_vel > th));
data_vel_cond(2).mus_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A39(sesh_num).tri_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th),...
    muscimol_alt_data.A39(sesh_num).bi_norm(muscimol_alt_data.A39(sesh_num).enc_vel > th));
data_vel_cond(3).mus_tribi_corr=corr_val(1,2);

sesh_num=17;
corr_val=corrcoef(muscimol_alt_data.A40(sesh_num).tri_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th),...
    muscimol_alt_data.A40(sesh_num).bi_norm(muscimol_alt_data.A40(sesh_num).enc_vel > th));
data_vel_cond(4).mus_tribi_corr=corr_val(1,2);

corr_val=corrcoef(muscimol_alt_data.A41(sesh_num).tri_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th),...
    muscimol_alt_data.A41(sesh_num).bi_norm(muscimol_alt_data.A41(sesh_num).enc_vel > th));
data_vel_cond(5).mus_tribi_corr=corr_val(1,2);



% w/ normalization to baseline

% plot fraction of session the mice spend moving wheel
% normed to baseline
mus_norm=vertcat(data_vel_cond(1:5).mus_perc)./vertcat(data_vel_cond(1:5).baseline_perc);
sal_norm=vertcat(data_vel_cond(1:5).saline_perc)./vertcat(data_vel_cond(1:5).baseline_perc);
figure; 
subplot(2,2,1); hold on;
plot(2,mean(mus_norm),'.k');
errorbar(2,mean(mus_norm),std(mus_norm)/sqrt(5),'k');
plot(1,mean(sal_norm),'.k');
errorbar(1,mean(sal_norm),std(sal_norm)/sqrt(5),'k');
plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Frac of session','wheel moving forward'});
title('Wheel engagement');
%ttest
[h,p]=ttest(mus_norm,sal_norm,'Tail','left'); %mus v sal
text(1.25,0.1,num2str(round(p,4)));
% example: ttest(x,y,'Tail','left');
% with 1-tailed, you're testing that the mean of x is LESS THAN the mean of y
%p=0.0212


% find mean wheel velocity during wheel movement
% normed to baseline
mus_norm=vertcat(data_vel_cond(1:5).mus_mean)./vertcat(data_vel_cond(1:5).baseline_mean);
sal_norm=vertcat(data_vel_cond(1:5).saline_mean)./vertcat(data_vel_cond(1:5).baseline_mean);

subplot(2,2,2); hold on;
plot(2,mean(mus_norm),'.k');
errorbar(2,mean(mus_norm),std(mus_norm)/sqrt(5),'k');
plot(1,mean(sal_norm),'.k');
errorbar(1,mean(sal_norm),std(sal_norm)/sqrt(5),'k');

plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Mean wheel velocity','when engaged (cm/s)'});
title('Velocity when wheel is engaged');
%ttest
[h,p]=ttest(mus_norm,sal_norm,'Tail','left'); %mus v sal
text(1.25,0.1,num2str(round(p,4)));
%p=0.0492

offset=0.1;
% normed to baseline
mus_tri_norm=vertcat(data_vel_cond(1:5).mus_mean_tri)./vertcat(data_vel_cond(1:5).base_mean_tri);
sal_tri_norm=vertcat(data_vel_cond(1:5).sal_mean_tri)./vertcat(data_vel_cond(1:5).base_mean_tri);
mus_bi_norm=vertcat(data_vel_cond(1:5).mus_mean_bi)./vertcat(data_vel_cond(1:5).base_mean_bi);
sal_bi_norm=vertcat(data_vel_cond(1:5).sal_mean_bi)./vertcat(data_vel_cond(1:5).base_mean_bi);

subplot(2,2,3); hold on;
plot(2-offset,mean(mus_tri_norm),'.k');
errorbar(2-offset,mean(mus_tri_norm),std(mus_tri_norm)/sqrt(5),'k');
plot(1-offset,mean(sal_tri_norm),'.k');
errorbar(1-offset,mean(sal_tri_norm),std(sal_tri_norm)/sqrt(5),'k');
plot(2+offset,mean(mus_bi_norm),'.m');
errorbar(2+offset,mean(mus_bi_norm),std(mus_bi_norm)/sqrt(5),'m');
plot(1+offset,mean(sal_bi_norm),'.m');
errorbar(1+offset,mean(sal_bi_norm),std(sal_bi_norm)/sqrt(5),'m');
plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Tri (k) & bi (m) activity','(norm) when muscles engaged'});
%ttest
[h,p]=ttest(mus_tri_norm,sal_tri_norm,'Tail','left'); %mus v sal
text(1,0.4,strcat('tri: ',num2str(round(p,4))));

[h,p]=ttest(mus_bi_norm,sal_bi_norm,'Tail','left'); %mus v sal
text(1,0.2,strcat('bi: ',num2str(round(p,4))));



%% cocon only when muscles active (>10% baseline max)

th=0.1;

% find active EMG periods and take CCI
%find times when tri or bi is above thresh, then take corresponding CCI
data_active=struct;
sesh_dur=1200000;
quant_val=0.9999;

% add baselines sessions 1-3, then average across them

% baseline 1
sesh_num=13;
mouse='A37';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(1).base1_mean=mean(vec); %mean CCI
data_active_bases(1).base1_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(1).base1_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(1).base1_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(1).base1_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(1).base1_tribi_corr=corr_vals(1,2);

mouse='A38';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(2).base1_mean=mean(vec); %mean CCI
data_active_bases(2).base1_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(2).base1_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(2).base1_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(2).base1_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(2).base1_tribi_corr=corr_vals(1,2);

mouse='A39';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(3).base1_mean=mean(vec); %mean CCI
data_active_bases(3).base1_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(3).base1_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(3).base1_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(3).base1_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(3).base1_tribi_corr=corr_vals(1,2);

mouse='A40';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(4).base1_mean=mean(vec); %mean CCI
data_active_bases(4).base1_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(4).base1_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(4).base1_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(4).base1_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(4).base1_tribi_corr=corr_vals(1,2);

mouse='A41';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(5).base1_mean=mean(vec); %mean CCI
data_active_bases(5).base1_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(5).base1_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(5).base1_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(5).base1_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(5).base1_tribi_corr=corr_vals(1,2);

% baseline 2
sesh_num=14;
mouse='A37';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(1).base2_mean=mean(vec); %mean CCI
data_active_bases(1).base2_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(1).base2_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(1).base2_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(1).base2_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(1).base2_tribi_corr=corr_vals(1,2);

mouse='A38';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(2).base2_mean=mean(vec); %mean CCI
data_active_bases(2).base2_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(2).base2_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(2).base2_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(2).base2_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(2).base2_tribi_corr=corr_vals(1,2);

mouse='A39';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(3).base2_mean=mean(vec); %mean CCI
data_active_bases(3).base2_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(3).base2_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(3).base2_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(3).base2_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(3).base2_tribi_corr=corr_vals(1,2);

mouse='A40';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(4).base2_mean=mean(vec); %mean CCI
data_active_bases(4).base2_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(4).base2_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(4).base2_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(4).base2_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(4).base2_tribi_corr=corr_vals(1,2);

mouse='A41';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(5).base2_mean=mean(vec); %mean CCI
data_active_bases(5).base2_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(5).base2_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(5).base2_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(5).base2_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(5).base2_tribi_corr=corr_vals(1,2);

% baseline 3
sesh_num=15;
mouse='A37';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(1).base3_mean=mean(vec); %mean CCI
data_active_bases(1).base3_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(1).base3_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(1).base3_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(1).base3_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(1).base3_tribi_corr=corr_vals(1,2);

mouse='A38';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(2).base3_mean=mean(vec); %mean CCI
data_active_bases(2).base3_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(2).base3_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(2).base3_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(2).base3_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(2).base3_tribi_corr=corr_vals(1,2);

mouse='A39';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(3).base3_mean=mean(vec); %mean CCI
data_active_bases(3).base3_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(3).base3_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(3).base3_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(3).base3_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(3).base3_tribi_corr=corr_vals(1,2);

mouse='A40';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(4).base3_mean=mean(vec); %mean CCI
data_active_bases(4).base3_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(4).base3_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(4).base3_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(4).base3_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(4).base3_tribi_corr=corr_vals(1,2);

mouse='A41';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active_bases(5).base3_mean=mean(vec); %mean CCI
data_active_bases(5).base3_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
data_active_bases(5).base3_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active_bases(5).base3_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active_bases(5).base3_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active_bases(5).base3_tribi_corr=corr_vals(1,2);

% average the baseline values together
for i=1:5
data_active(i).baseline_mean=mean([data_active_bases(i).base1_mean data_active_bases(i).base2_mean data_active_bases(i).base3_mean]); %mean CCI
data_active(i).baseline_perc=mean([data_active_bases(i).base1_perc data_active_bases(i).base2_perc data_active_bases(i).base3_perc]);
data_active(i).baseline_tribi_pdist=mean([data_active_bases(i).base1_tribi_pdist data_active_bases(i).base2_tribi_pdist data_active_bases(i).base3_tribi_pdist]);
data_active(i).baseline_tri_mean=mean([data_active_bases(i).base1_tri_mean data_active_bases(i).base2_tri_mean data_active_bases(i).base3_tri_mean]);
data_active(i).baseline_bi_mean=mean([data_active_bases(i).base1_bi_mean data_active_bases(i).base2_bi_mean data_active_bases(i).base3_bi_mean]);
data_active(i).baseline_tribi_corr=mean([data_active_bases(i).base1_tribi_corr data_active_bases(i).base2_tribi_corr data_active_bases(i).base3_tribi_corr]);
end


%muscimol
sesh_num=16;
mouse='A37';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(1).mus=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(1).mus_mean=mean(vec);
data_active(1).mus_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(1).mus_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(1).mus_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(1).mus_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(1).mus_tribi_corr=corr_vals(1,2);

mouse='A38';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(2).mus=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
data_active(2).mus_mean=mean(muscimol_cc_data.(mouse)(sesh_num).CCI(inds));
data_active(2).mus_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(2).mus_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(2).mus_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(2).mus_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(2).mus_tribi_corr=corr_vals(1,2);

mouse='A39';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(3).mus=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(3).mus_mean=mean(vec);
data_active(3).mus_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(3).mus_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(3).mus_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(3).mus_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(3).mus_tribi_corr=corr_vals(1,2);

sesh_num=18;
mouse='A40';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(4).mus=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
data_active(4).mus_mean=mean(muscimol_cc_data.(mouse)(sesh_num).CCI(inds));
data_active(4).mus_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(4).mus_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(4).mus_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(4).mus_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(4).mus_tribi_corr=corr_vals(1,2);

mouse='A41';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(5).mus=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(5).mus_mean=mean(vec);
data_active(5).mus_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(5).mus_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(5).mus_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(5).mus_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(5).mus_tribi_corr=corr_vals(1,2);

%saline
sesh_num=18;
mouse='A37';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(1).sal=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(1).sal_mean=mean(vec);
data_active(1).sal_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(1).sal_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(1).sal_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(1).sal_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(1).sal_tribi_corr=corr_vals(1,2);

mouse='A38';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(2).sal=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
data_active(2).sal_mean=mean(muscimol_cc_data.(mouse)(sesh_num).CCI(inds));
data_active(2).sal_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(2).sal_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(2).sal_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(2).sal_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(2).sal_tribi_corr=corr_vals(1,2);

mouse='A39';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(3).sal=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(3).sal_mean=mean(vec);
data_active(3).sal_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(3).sal_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(3).sal_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(3).sal_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(3).sal_tribi_corr=corr_vals(1,2);

sesh_num=16;
mouse='A40';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(4).sal=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
data_active(4).sal_mean=mean(muscimol_cc_data.(mouse)(sesh_num).CCI(inds));
data_active(4).sal_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(4).sal_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(4).sal_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(4).sal_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(4).sal_tribi_corr=corr_vals(1,2);

mouse='A41';
inds=logical((muscimol_cc_data.(mouse)(sesh_num).tri_norm>th)+(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
data_active(5).sal=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec=muscimol_cc_data.(mouse)(sesh_num).CCI(inds);
vec(isnan(vec))=[];
data_active(5).sal_mean=mean(vec);
data_active(5).sal_perc=length(muscimol_cc_data.(mouse)(sesh_num).CCI(inds))/sesh_dur;
active_tri=muscimol_cc_data.(mouse)(sesh_num).tri_norm(inds);
active_bi=muscimol_cc_data.(mouse)(sesh_num).bi_norm(inds);
data_active(5).sal_tribi_pdist=pdist([active_tri; active_bi])/length(active_tri); %tri-bi pdist when muscles active
data_active(5).sal_tri_mean=mean(muscimol_cc_data.(mouse)(sesh_num).tri_norm(muscimol_cc_data.(mouse)(sesh_num).tri_norm>th));
data_active(5).sal_bi_mean=mean(muscimol_cc_data.(mouse)(sesh_num).bi_norm(muscimol_cc_data.(mouse)(sesh_num).bi_norm>th));
corr_vals=corrcoef(active_tri,active_bi);
data_active(5).sal_tribi_corr=corr_vals(1,2);

% normed to baseline
mus_norm=vertcat(data_active(1:5).mus_perc)./vertcat(data_active(1:5).baseline_perc);
sal_norm=vertcat(data_active(1:5).sal_perc)./vertcat(data_active(1:5).baseline_perc);

figure; 
subplot(2,2,1); hold on;
plot(2,mean(mus_norm),'.k');
errorbar(2,mean(mus_norm),std(mus_norm)/sqrt(5),'k');
plot(1,mean(sal_norm),'.k');
errorbar(1,mean(sal_norm),std(sal_norm)/sqrt(5),'k');
plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Fraction of session','with active tri or bi','normed to baseline'});
%ttest
[h,p]=ttest(mus_norm,sal_norm,'Tail','left'); %mus v sal
text(1.25,0.2,num2str(round(p,4)));


% normed to baseline
mus_norm=vertcat(data_active(1:5).mus_mean)./vertcat(data_active(1:5).baseline_mean);
sal_norm=vertcat(data_active(1:5).sal_mean)./vertcat(data_active(1:5).baseline_mean);

subplot(2,2,2); hold on;
plot(2,mean(mus_norm),'.k');
errorbar(2,mean(mus_norm),std(mus_norm)/sqrt(5),'k');
plot(1,mean(sal_norm),'.k');
errorbar(1,mean(sal_norm),std(sal_norm)/sqrt(5),'k');
plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Mean CCI over active','muscle periods'});
%ttest
[h,p]=ttest(mus_norm,sal_norm,'Tail','left'); %mus v sal
text(1.25,0.2,num2str(round(p,4)));

% normed to baseline
mus_tri_norm=vertcat(data_active(1:5).mus_tri_mean)./vertcat(data_active(1:5).baseline_tri_mean);
sal_tri_norm=vertcat(data_active(1:5).sal_tri_mean)./vertcat(data_active(1:5).baseline_tri_mean);
mus_bi_norm=vertcat(data_active(1:5).mus_bi_mean)./vertcat(data_active(1:5).baseline_bi_mean);
sal_bi_norm=vertcat(data_active(1:5).sal_bi_mean)./vertcat(data_active(1:5).baseline_bi_mean);

subplot(2,2,3); hold on;
plot(2-offset,mean(mus_tri_norm),'.k');
errorbar(2-offset,mean(mus_tri_norm),std(mus_tri_norm)/sqrt(5),'k');
plot(1-offset,mean(sal_tri_norm),'.k');
errorbar(1-offset,mean(sal_tri_norm),std(sal_tri_norm)/sqrt(5),'k');
plot(2+offset,mean(mus_bi_norm),'.m');
errorbar(2+offset,mean(mus_bi_norm),std(mus_bi_norm)/sqrt(5),'m');
plot(1+offset,mean(sal_bi_norm),'.m');
errorbar(1+offset,mean(sal_bi_norm),std(sal_bi_norm)/sqrt(5),'m');
plot([0 2.5],[1 1],'--k');
xlim([0.5 2.5]);
ylim([0 1.5]);
box off;
axis square;
xticks(1:2);
xticklabels({'saline','muscimol'});
ylabel({'Tri (k) & bi (m) activity','(norm) when muscles engaged'});
%ttest
[h,p]=ttest(mus_tri_norm,sal_tri_norm,'Tail','left'); %mus v sal
text(1,0.4,strcat('tri: ',num2str(round(p,4))));
[h,p]=ttest(mus_bi_norm,sal_bi_norm,'Tail','left'); %mus v sal
text(1,0.2,strcat('bi: ',num2str(round(p,4))));

suptitle('During times tri OR bi above quiescence threshold (10% max)');
