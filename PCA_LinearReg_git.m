
%% PCA and Linear Regression

% Run this matlab script when the folder specified in data_dir_root is open, which should
% contain the processed data

for mouse = 1:4
    
    %% Data Directory
    data_dir_root='G:\manuscript_data\';  
    
    %% Which Mouse?
    mouse
    if mouse==1
        mouse_name='A32';
    elseif mouse==2
        mouse_name='A33';
    elseif mouse==3
        mouse_name='A35';
    elseif mouse==4
        mouse_name='A36';
    end
    
    saveon=1;
%     load_lambda=1;
%     save_lambda=0;
%     load_lambda_joint=1;
%     save_lambda_joint=0;
%     load_lambda_joint_scale=0;
%     save_lambda_joint_scale=1;
%     cross_validate=1;

    
    %% What to plot?
    plot_lin_reg_traces=0;
    plot_lin_reg_traces_joint=0;
    plot_lin_reg_weights_heatmap=0;
    plot_in_pyr_figs=0;
    plot_cleaned_neuron_examples=0;
    
    %% Set Parameters
    downsample_int=30;
    trial_idxs=101:400;
    lw=1.5; %this is line width for plots
    time_shift=0;
    
    %lin reg params
    numval_trials=10; %MUST BE 10 in ver 15 and after
    num_perms_lambda=500; %(must be at least 2) %num iterations of testing for R2 (on diff permutations)
    
    %exclusion criteria
    minfr=1; % 3, Hz (Set to 0 to turn off, ie for making cum sum of FRs fro fig 2)
    slopeThresh=0.2; %0.3, MAX value for slope (Set to inf to turn off) This is the pdist change on the waveform
    fracfrFactor=1.5; %(Set to inf to turn off)
    
    %PCA criteria
    numPCs=6; %for PCA %had been 4 (tried 6 for comparison with EMG)
    numPCs_emg=6;
    
    %similarity index alignment
    PCAdim = 6; %had been 4 (tried 6 for comparison with EMG)
    
    %cleaning criterion
    numPCs_forClean=6;
    
    %subtypes criteria
    width_IN=0.4167; %UPDATED 6/10/19 based on vgat data
    width_pyr=0.4167;
    
    if strcmp('A32',mouse_name)
        drop_sessions=[4,9];  % new hard th exclusions
    elseif strcmp('A33',mouse_name)
        drop_sessions=[1,2,4,5,6,8,10,11,16,17,18];
    elseif strcmp('A35',mouse_name)
        drop_sessions=[2:8,15,17]; % new hard th exclusions
    elseif strcmp('A36',mouse_name)
        drop_sessions=[4,7,9,10]; % new hard th exclusions
    end
    
    
    %% Load Data and Parameters
    samplingrate=30000;
    ms2samp=0.001*samplingrate/downsample_int;
    sr=0.001*ms2samp;
    
    load(strcat(data_dir_root,mouse_name,'_emg_spikes.mat'));
    
    load(strcat(data_dir_root,mouse_name,'_widths_tot_1point5ms_FILTERED.mat'));
    load(strcat(data_dir_root,mouse_name,'_allsesh_waveform_data_filtered.mat'));
    load(strcat(data_dir_root,mouse_name,'_waveform_stability_data'));
    
    waveform_data=eval(strcat(mouse_name,'_allsesh_waveform_data_filtered'));
    slopes_data=data;
    
    % correct errors
    % in A35 session 6, bug where only 2 units  of 23 load in save spike times,
    % for now, delete the 21 that are not showing up,
    % ie keep # 1 (64) and # 7 (193)
    if strcmp('A35',mouse_name)
        waveform_data(6).meanWF=[waveform_data(6).meanWF(1,:); waveform_data(6).meanWF(7,:)];
        widths_new([152:156 158:173])=[];
    end
    
    %patch uneven number
    %in A36 session 3, unit 264 is not included in the waveform and width calc
    if strcmp('A36',mouse_name)
        mouse_data(3).spikes(12)=[];
        mouse_data(3).firingrate(12)=[];
    end
    %note, A36, session 3, there is 1 waveform missing (total neurons=34 but
    %total waveforms=33)
    
    if strcmp('A33',mouse_name)
        width_15ms=widths_new;
    elseif strcmp('A35',mouse_name)
        width_15ms=widths_new;
    end
    
    %% drop bad sessions
    
    %create vector with session number and neuron number of each neuron
    neuron_session_ID=[];
    neuron_number_ID=[];
    for i = 1:length(mouse_data)
        neuron_session_ID=[neuron_session_ID; i*ones(length(mouse_data(i).spikes),1)];
        neuron_number_ID=[neuron_number_ID 1:length(mouse_data(i).spikes)];
    end
    neuron_origin_ID=[neuron_session_ID neuron_number_ID'];
    
    if ~isempty(drop_sessions)
        mouse_data(drop_sessions)=[];
        waveform_data(drop_sessions)=[];
        slopes_data(drop_sessions)=[];
        throw_vec=zeros(length(neuron_session_ID),1);
        for i=1:length(drop_sessions)
            throw_vec=throw_vec+(neuron_session_ID==drop_sessions(i));
        end
        throw_vec=logical(throw_vec);
        width_15ms(throw_vec)=[];
        neuron_session_ID(throw_vec)=[];
        neuron_number_ID(throw_vec)=[];
        neuron_origin_ID(throw_vec,:)=[];
    end
    
    %% initialization
    
    %list conditions and muscles?
    conds={'cc','alt1','alt2','inact'};
    muscle_names={'bi','edc','pal','pec','spdel','triM'};
    
    %Add col with total neuron number
    neurons_total=zeros(1,length(mouse_data));
    for i=1:length(mouse_data)
        mouse_data(i).neurons_total=length(mouse_data(i).firingrate);
        neurons_total(i)=length(mouse_data(i).firingrate);
    end
    
    %Create Data Matrices
    numSessions=length(mouse_data);
    numMuscles=length(muscle_names);
    
    %% Implement Neuron Selection via FR/slopes
    
    dropNeurs=cell(numSessions,1);
    keepVecTot=[];
    neurons_select=zeros(1,length(mouse_data));
    totalfrCell=cell(numSessions,1);
    
    for j=1:numSessions
        keepVec=ones(1,length(mouse_data(j).firingrate));
        totalfr=nan(length(conds)-1,length(mouse_data(j).firingrate));
        for i=1:3 %Only look at reg conds (not special extras)
            cond_trials=strcat('trials_',conds{i});
            for n=1:length(mouse_data(j).firingrate)
                totalfr(i,n)=(ms2samp/0.001)*sum(sum(mouse_data(j).spikes(n).(cond_trials)(:,trial_idxs)))/numel(mouse_data(j).spikes(n).(cond_trials)(:,trial_idxs));
            end
        end
        totalfrCell{j}=totalfr;
        
        %exclude neurons with FR<min in all conditions
        totalfr_summary=sum(totalfr>=minfr,1);
        dropNeurs{j}=find(totalfr_summary==0);
        for n=1:length(mouse_data(j).firingrate)
            if any(strcmp(conds,'alt1')) && any(strcmp(conds,'alt2'))
                if fracfrFactor<inf % factor change
                    if totalfr(strcmp(conds,'alt2'),n)>totalfr(strcmp(conds,'alt1'),n)*fracfrFactor || ...
                            totalfr(strcmp(conds,'alt2'),n)<totalfr(strcmp(conds,'alt1'),n)/fracfrFactor
                        dropNeurs{j}=[dropNeurs{j} n];
                    end
                end
            end
        end
        
        %use slope of pdist (WF stability) to exclude Ns
        dropNeurs{j}=[dropNeurs{j} find(slopes_data(j).slopes>=slopeThresh);];
        dropNeurs{j}=unique(dropNeurs{j});
        keepVec(dropNeurs{j})=0;
        keepVecTot=[keepVecTot keepVec];
        fprintf('%d neurons dropped in Session %d\n',length(dropNeurs{j}),j);
        
        mouse_data(j).firingrate_select=mouse_data(j).firingrate;
        mouse_data(j).spikes_select=mouse_data(j).spikes;
        mouse_data(j).firingrate_select(dropNeurs{j})=[];
        mouse_data(j).spikes_select(dropNeurs{j})=[];
        mouse_data(j).neurons_select=length(mouse_data(j).firingrate_select);
        neurons_select(j)=length(mouse_data(j).firingrate_select);
        mouse_data(j).waveform=waveform_data(j).meanWF;
        mouse_data(j).waveform_select=waveform_data(j).meanWF;
        mouse_data(j).waveform_select(dropNeurs{j},:)=[];
    end
    
    numNeurons=zeros(numSessions,1);
    numTrials=zeros(length(conds),numSessions);
    neuron_origin_ID_select=neuron_origin_ID(logical(keepVecTot)',:);
    
    % %make pie chart of selected neurons
    % if sum(neurons_total)-sum(neurons_select) >0
    %     figure;
    %     pie([sum(neurons_total)-sum(neurons_select) sum(neurons_select)],{strcat('Removed neurons:',' ',num2str((sum(neurons_total)-sum(neurons_select))/sum(neurons_total))),...
    %         strcat('Selected neurons:',' ',num2str(sum(neurons_select)/sum(neurons_total)))});
    % end
    
    %make mat of all waveforms
    waveforms_select=[];
    for i=1:numSessions
        waveforms_select=[waveforms_select; mouse_data(i).waveform_select];
    end
    
    %How many neurons and trials
    sess2neurs=cell(numSessions,1);
    for j=1:numSessions
        numNeurons(j)=length(mouse_data(j).firingrate_select);
        for i=1:length(conds)
            cond_trials=strcat('trials_',conds{i});
            % Spike data(i) is numTrials x numSamples (T)
            for n=1:numNeurons(j)
                numTrials(i,j)=size(mouse_data(j).firingrate_select(n).(cond_trials)(:,trial_idxs),1);
            end
        end
        sess2neurs{j}=sum(numNeurons(1:j-1))+(1:numNeurons(j));
    end
    totneurs=sum(numNeurons(:));
    fprintf('\n%d total neurons in all sessions\n',totneurs);

    %% Build Data matrices
    Data=struct;
    perm_trials=nan(num_perms_lambda,numTrials(1,1));
    if numval_trials==1
        % for leave-one-out
        for i=1:numTrials(1,1)
            if i==1
                perm_trials(i,:)=i:numTrials(1,1);
            else
                perm_trials(i,:)=[i:numTrials(1,1) 1:i-1];
            end
        end
    else
        % Do permutations many times
        for i=1:num_perms_lambda
            perm_trials(i,:)=randperm(numTrials(1,1));
        end
    end
    
    for i=1:length(conds)
        
        %initialize variables
        cond_trials=strcat('trials_',conds{i});
        mean_fr=nan(length(trial_idxs),totneurs); % size TxN
        mean_emg=zeros(length(trial_idxs),length(muscle_names)); % size TxN
        
        train_trials=cell(numSessions,num_perms_lambda);
        val_trials=cell(numSessions,num_perms_lambda);
        train_neur=cell(num_perms_lambda,1);
        val_neur=cell(num_perms_lambda,1);
        train_emg=cell(num_perms_lambda,1);
        val_emg=cell(num_perms_lambda,1);
        % new
        test_trials=cell(numSessions,num_perms_lambda);
        test_neur=cell(num_perms_lambda,1);
        test_emg=cell(num_perms_lambda,1);
        
        %makes trial averages of neural and emg (for use after validation)
        for j=1:numSessions
            for n=1:numNeurons(j)
                mean_fr(:,sess2neurs{j}(n))=1000*mean(mouse_data(j).firingrate_select(n).(cond_trials)(:,trial_idxs-time_shift),1)'; % in Hz
            end
            for m=1:length(muscle_names)
                mean_emg(:,m)=mean_emg(:,m)+sum(mouse_data(j).emg.(cond_trials).(muscle_names{m})(:,trial_idxs),1)';
            end
        end
        mean_emg=mean_emg/sum(numTrials(i,:));
        
        for b=1:num_perms_lambda
            train_neur{b}=nan(length(trial_idxs-time_shift),totneurs);
            train_emg{b}=zeros(length(trial_idxs),length(muscle_names));
            val_neur{b}=nan(length(trial_idxs-time_shift),totneurs);
            val_emg{b}=zeros(length(trial_idxs),length(muscle_names));
            % new
            test_neur{b}=nan(length(trial_idxs-time_shift),totneurs);
            test_emg{b}=zeros(length(trial_idxs),length(muscle_names));
            
            for j=1:numSessions
                val_trials{j,b}=perm_trials(b,1:10);
                train_trials{j,b}=perm_trials(b,11:20);
                test_trials{j,b}=perm_trials(b,21:30);
                
                for n=1:numNeurons(j)
                    train_neur{b}(:,sess2neurs{j}(n))=1000*mean(mouse_data(j).firingrate_select(n).(cond_trials)(train_trials{j,b},trial_idxs-time_shift),1)';
                    val_neur{b}(:,sess2neurs{j}(n))=1000*mean(mouse_data(j).firingrate_select(n).(cond_trials)(val_trials{j,b},trial_idxs-time_shift),1)';
                    test_neur{b}(:,sess2neurs{j}(n))=1000*mean(mouse_data(j).firingrate_select(n).(cond_trials)(test_trials{j,b},trial_idxs-time_shift),1)';
                end
                for m=1:length(muscle_names)
                    train_emg{b}(:,m)=train_emg{b}(:,m)+sum(mouse_data(j).emg.(cond_trials).(muscle_names{m})(train_trials{j,b},trial_idxs),1)';
                    val_emg{b}(:,m)=val_emg{b}(:,m)+sum(mouse_data(j).emg.(cond_trials).(muscle_names{m})(val_trials{j,b},trial_idxs),1)';
                    test_emg{b}(:,m)=test_emg{b}(:,m)+sum(mouse_data(j).emg.(cond_trials).(muscle_names{m})(test_trials{j,b},trial_idxs),1)';
                end
            end
            train_emg{b}=train_emg{b}/sum(length(horzcat(train_trials{:,b})));
            val_emg{b}=val_emg{b}/sum(length(horzcat(val_trials{:,b})));
            test_emg{b}=test_emg{b}/sum(length(horzcat(test_trials{:,b})));
        end
        
        Data(i).fr=mean_fr;
        Data(i).emg=mean_emg;
        Data(i).train_trials=train_trials;
        Data(i).val_trials=val_trials;
        Data(i).train_neur=train_neur;
        Data(i).train_emg=train_emg;
        Data(i).val_neur=val_neur;
        Data(i).val_emg=val_emg;
        Data(i).times=0.001*ms2samp*(1:length(trial_idxs));
        times = sr*(0:length(trial_idxs)-1);
        %new
        Data(i).test_trials=test_trials;
        Data(i).test_neur=test_neur;
        Data(i).test_emg=test_emg;
    end
    
    % figure;
    % for i=1:3, plot(mean(Data(i).fr,2),'linewidth',2); hold on; end
    % set(gca,'fontsize',16); legend(conds);
    % xlabel('time'); ylabel('Mean Firing Rate');
    
   

    %% Look at subtype FRs

    width_15ms_keep=width_15ms(logical(keepVecTot)');
    
    Data_subtypes=struct;
    Data_subtypes(1).condition = 'cc';
    Data_subtypes(2).condition = 'alt1';
    Data_subtypes(3).condition = 'alt2';
    Data_subtypes(4).condition = 'inact';
    Data_subtypes(5).condition = 'waveform';
    for i=1:length(conds)
        Data_subtypes(i).interneurons = Data(i).fr(:,width_15ms_keep<width_IN);
        Data_subtypes(i).pyramidals = Data(i).fr(:,width_15ms_keep>width_pyr);
    end
    Data_subtypes(5).interneurons = waveforms_select(width_15ms_keep<width_IN,:)';
    Data_subtypes(5).pyramidals = waveforms_select(width_15ms_keep>width_pyr,:)';
    
    if plot_in_pyr_figs %make IN vs Pyr plots
        input=Data_subtypes;
        xvec=1:300;
        xlim_vals=[0 301];
        if strcmp(mouse_name,'A32')
            ylim_vals1=[0 80];
            ylim_vals2=[0 20];
        elseif strcmp(mouse_name,'A33')
            ylim_vals1=[0 90];
            ylim_vals2=[0 25];
        elseif strcmp(mouse_name,'A35')
            ylim_vals1=[0 90];
            ylim_vals2=[0 25];
        elseif strcmp(mouse_name,'A36')
            ylim_vals1=[0 90];
            ylim_vals2=[0 25];
        end
        
        figure;
        subplot(2,4,1);
        shadedErrorBar(xvec,mean(input(4).interneurons'),std(input(4).interneurons')/sqrt(size(input(4).interneurons',1)),{'Color',[1 0 0],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals1);
        box off;
        ylabel('Firing rate (spks/s)');
        title('Inact');
        subplot(2,4,2);
        shadedErrorBar(xvec,mean(input(2).interneurons'),std(input(2).interneurons')/sqrt(size(input(2).interneurons',1)),{'Color',[1 0 0],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals1);
        box off;
        title('Alt1');
        subplot(2,4,3);
        shadedErrorBar(xvec,mean(input(1).interneurons'),std(input(1).interneurons')/sqrt(size(input(1).interneurons',1)),{'Color',[1 0 0],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals1);
        box off;
        title('Cocon');
        subplot(2,4,4);
        shadedErrorBar(xvec,mean(input(3).interneurons'),std(input(3).interneurons')/sqrt(size(input(3).interneurons',1)),{'Color',[1 0 0],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals1);
        box off;
        title('Alt2');
        
        subplot(2,4,5);
        shadedErrorBar(xvec,mean(input(4).pyramidals'),std(input(4).pyramidals')/sqrt(size(input(4).pyramidals',1)),{'Color',[0 0 1],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals2);
        box off;
        ylabel('Firing rate (spks/s)');
        subplot(2,4,6);
        shadedErrorBar(xvec,mean(input(2).pyramidals'),std(input(2).pyramidals')/sqrt(size(input(2).pyramidals',1)),{'Color',[0 0 1],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals2);
        box off;
        subplot(2,4,7);
        shadedErrorBar(xvec,mean(input(1).pyramidals'),std(input(1).pyramidals')/sqrt(size(input(1).pyramidals',1)),{'Color',[0 0 1],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals2);
        box off;
        subplot(2,4,8);
        shadedErrorBar(xvec,mean(input(3).pyramidals'),std(input(3).pyramidals')/sqrt(size(input(3).pyramidals',1)),{'Color',[0 0 1],'LineWidth',2},1);
        xlim(xlim_vals);
        ylim(ylim_vals2);
        box off;
        xlabel('Time (ms)');
        
        titleStr = '%d INs and %d PYRs';
        suptitle(sprintf(titleStr,size(input(4).interneurons',1),size(input(4).pyramidals',1)));
        
        % make normalized plots
        normed_IN=(mean(input(1).interneurons')-min(mean(input(1).interneurons')))/(max(mean(input(1).interneurons'))-min(mean(input(1).interneurons')));
        normed_Pyr=(mean(input(1).pyramidals')-min(mean(input(1).pyramidals')))/(max(mean(input(1).pyramidals'))-min(mean(input(1).pyramidals')));
        figure; hold on;
        plot(xvec,normed_IN,'r');
        plot(xvec,normed_Pyr,'b');
        xlim(xlim_vals);
        title(strcat(mouse_name,'-Cocon normed IN (red) and PYR (blue)'));
    end
    
    %% break Data into subtypes for lin reg
    Data_pyr=Data;
%     Data_in=Data;
    for i=1:length(conds)
        Data_pyr(i).fr(:,width_15ms_keep<width_pyr)=[];
        for ii=1:length(Data(i).train_neur)
            Data_pyr(i).train_neur{ii}(:,width_15ms_keep<width_pyr)=[];
            Data_pyr(i).val_neur{ii}(:,width_15ms_keep<width_pyr)=[];
            Data_pyr(i).test_neur{ii}(:,width_15ms_keep<width_pyr)=[];
        end
%             Data_in(i).fr(:,width_15ms_keep>width_IN)=[];
%             for ii=1:length(Data(i).train_neur)
%                 Data_in(i).train_neur{ii}(:,width_15ms_keep>width_IN)=[];
%                 Data_in(i).val_neur{ii}(:,width_15ms_keep>width_IN)=[];
%             end
    end
    


 %% try to take PCs of alt then orthogonalize all data of cocon w/ respect to those
 % PCs
 
 %all neurons
 input=Data;
 fr_mat=input(2).fr; %alt1
 % mean center, only subtract off total mean for each neuron across the conditions
 fr_mat=fr_mat-repmat(mean(fr_mat),size(fr_mat,1),1);
 [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
 PCs=PC;
 % get subspace projection of top 6 PCs
 projs=PCs(:,1:6)*PCs(:,1:6)';
 % multiply by identity matrix to get orthogonal complement
 ortho_com = eye(size(projs,1))*projs;
 %find the cocon activity that is orthogonal to the alt subspace
 cc_act_ortho = ortho_com * input(1).fr';

 %pyr neurons
 input=Data_pyr;
 fr_mat=input(2).fr; %alt1
 % mean center, only subtract off total mean for each neuron across the conditions
 fr_mat=fr_mat-repmat(mean(fr_mat),size(fr_mat,1),1);
 [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
 PCs=PC;
 % get subspace projection of top 6 PCs
 projs=PCs(:,1:6)*PCs(:,1:6)';
 % multiply by identity matrix to get orthogonal complement
 ortho_com = eye(size(projs,1))*projs;
 %find the cocon activity that is orthogonal to the alt subspace
 cc_act_ortho_pyr = ortho_com * input(1).fr';

    %% Linear Regression on CC and A1 models
    
    lambda_vec=0.5:0.01:2;
    %lambda_vec=-4:0.02:2;
    warning('off','stats:regress:RankDefDesignMat');
    
    Test_R2_ortho_cc=nan*ones(num_perms_lambda,1);
    Test_R2_ortho_cc_pyr=nan*ones(num_perms_lambda,1);
    
    
     for k = 1:2 %1 for all neurons, 2 for just pyr
        if k == 1
            input=Data;
            disp('lin reg all');
        elseif k == 2
            input=Data_pyr;
            disp('lin reg pyr');
        end
        R2_conds=nan*ones(length(conds));
        Y_train_opt=cell(length(conds),1);
        Y=cell(length(conds),length(conds));
        Y_fit=cell(length(conds),length(conds));

        for i=1:length(conds)-1 %only do CC and Alt1, and ALSO A2 for the weights
            %if ~load_lambda
            %Train_R2=nan*ones(length(lambda_vec),num_perms_lambda);
            Train_R2=nan*ones(length(lambda_vec),1);
            %Val_R2=nan*ones(length(lambda_vec),num_perms_lambda);
            Val_R2=nan*ones(length(lambda_vec),1);
            M_local=cell(length(lambda_vec));
            % new
            Test_R2=nan*ones(num_perms_lambda,length(conds));
            Y_fit=cell(num_perms_lambda,length(conds));
            Y_real=cell(num_perms_lambda,length(conds));
            lambda_opt=nan*ones(1,num_perms_lambda);
            M=cell(num_perms_lambda);
            
            
            
            % here, lambda calculated from train and val and tested on test trials in EACH permutation
            for b=1:num_perms_lambda
                %b
                % Fit the model on the training set of that condition
                X_train=input(i).train_neur{b};
                X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
                Y_train=input(i).train_emg{b};
                Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
                X_val=input(i).val_neur{b};
                X_val=X_val-repmat(mean(X_val,1),size(X_val,1),1); % mean center X_val
                Y_val=input(i).val_emg{b};
                Y_val=Y_val-repmat(mean(Y_val,1),size(Y_val,1),1); % mean center Y_val
                
                for j=1:length(lambda_vec)
                    [Train_R2(j), M_local{j}, ~] = anyRegress( X_train, Y_train, 'ridge',10^lambda_vec(j));
                    [Val_R2(j),~,~] = anyRegress( X_val, Y_val, 'M',M_local{j} );
                end
                % pick best lambda
                [max_R2,lambda_opt_ind]=max(Val_R2); %find lambda that gave highest R2
                lambda_opt(b)=10^(lambda_vec(lambda_opt_ind));
                
                %now combine the test and val trials
                X_train=(input(i).train_neur{b}+input(i).val_neur{b})./2; %mean of test and val
                %X_train=[X_train ones(row_num,1)]; %add col of ones for offset instead of mean center
                X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
                Y_train=(input(i).train_emg{b}+input(i).val_emg{b})./2;%mean of test and val
                Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
                
                % get the M for train data, use within-iteration lambda
                [~, M{b}, ~] = anyRegress(X_train, Y_train, 'ridge',lambda_opt(b));
                
                % test the weights
                for c=1:length(conds)
                    % use only subset of trials (test)
                    X=input(c).test_neur{b};
                    X=X-repmat(mean(X,1),size(X,1),1); % mean center X_val
                    Y=[input(c).test_emg{b}];
                    Y=Y-repmat(mean(Y,1),size(Y,1),1); % mean center Y_train
                    Y_real{b,c}=Y; %record your real emg for examples
                    [Test_R2(b,c),~,Y_fit{b,c}]=anyRegress(X, Y, 'M',M{b} );
                end
                
                %test on alt1 model on ortho cocon
                if k==1 && i==2
                    % use only subset of trials (test)
                    X=cc_act_ortho';
                    X=X-repmat(mean(X,1),size(X,1),1); % mean center X_val
                    Y=[input(1).test_emg{b}]; %use same EMG as regular cc
                    Y=Y-repmat(mean(Y,1),size(Y,1),1); % mean center Y_train
                    [Test_R2_ortho_cc(b),~,~]=anyRegress(X, Y, 'M',M{b} );
                end
                
                if k==2 && i==2
                    % use only subset of trials (test)
                    X=cc_act_ortho_pyr';
                    X=X-repmat(mean(X,1),size(X,1),1); % mean center X_val
                    Y=[input(1).test_emg{b}]; %use same EMG as regular cc
                    Y=Y-repmat(mean(Y,1),size(Y,1),1); % mean center Y_train
                    [Test_R2_ortho_cc_pyr(b),~,~]=anyRegress(X, Y, 'M',M{b} );
                end
            end
            
            %average over permutations
            Test_R2_mean=mean(Test_R2); %average over permutations
            Test_R2_ortho_cc_mean=mean(Test_R2_ortho_cc);
            Test_R2_ortho_cc_pyr_mean=mean(Test_R2_ortho_cc_pyr);
            
            %identify the iteration whose R2 vals were closest in aggregate to mean
            Test_R2_diff=nan*ones(size(Test_R2));
            for j=1:size(Test_R2,1)
                Test_R2_diff(j,:)=abs(Test_R2(j,:)-Test_R2_mean);
            end
            Test_R2_diff_sum=sum(Test_R2_diff(:,1:3),2);
            [val,ind]=min(Test_R2_diff_sum);
            
            if i==1
                Y_real_CCmodel=Y_real(ind,:);
                Y_fit_CCmodel=Y_fit(ind,:);
                M_CCmodel=M;
                Test_R2_CCmodel=Test_R2_mean;
                lambda_opt_CCmodel=lambda_opt;
                
                Y_fit_CCmodel_Xperms=Y_fit;
            elseif i==2
                Y_real_A1model=Y_real(ind,:);
                Y_fit_A1model=Y_fit(ind,:);
                M_A1model=M;
                Test_R2_A1model=Test_R2_mean;
                lambda_opt_A1model=lambda_opt;
                
                Y_fit_A1model_Xperms=Y_fit;
            elseif i==3    
                M_A2model=M;
                
            end
        end

        
        figure;
        bar([Test_R2_CCmodel(1:3); Test_R2_A1model(1:3)]);
        set(gca,'fontsize',20)
        xlabel('Condition used for fitting M');
        ax=gca;
        set(ax,'XTickLabel',conds);
        legend(conds(1:3));
        if k == 1
            title('R^2 ALL neurons');
        elseif k == 2
            title('R^2 PYR neurons');
        end
        
        if plot_lin_reg_traces
            ylim_vals=[-0.2 0.26];
            for i=1:2 %only alt1 and cocon
                if i==1
                    Y_real_input=Y_real_CCmodel;
                    Y_fit_input=Y_fit_CCmodel;
                elseif i==2
                    Y_real_input=Y_real_A1model;
                    Y_fit_input=Y_fit_A1model;
                end
                
                figure;
                subplot(2,3,1); hold on;
                plot(Y_real_input{1,1}(:,1),'k');
                ylabel('bi');
                title('cc');
                plot(Y_fit_input{1,1}(:,1),'g');
                axis square;
                box off;
                
                subplot(2,3,2); hold on;
                plot(Y_real_input{1,2}(:,1),'k');
                title('alt1');
                plot(Y_fit_input{1,2}(:,1),'g');
                axis square;
                box off;

                subplot(2,3,3); hold on;
                plot(Y_real_input{1,3}(:,1),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,1),'g');
                axis square;
                box off;

                subplot(2,3,4); hold on;
                plot(Y_real_input{1,1}(:,6),'k');
                ylabel('tri');
                plot(Y_fit_input{1,1}(:,6),'g');
                axis square;
                box off;

                subplot(2,3,5); hold on;
                plot(Y_real_input{1,2}(:,6),'k');
                plot(Y_fit_input{1,2}(:,6),'g');
                axis square;
                box off;
                
                subplot(2,3,6); hold on;
                plot(Y_real_input{1,3}(:,6),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,6),'g');
                axis square;
                box off;
                
                if k == 1
                    if i==1
                        suptitle('Train on CC ALL');
                    elseif i==2
                        suptitle('Train on Alt1 ALL');
                    end
                elseif k == 2
                    if i==1
                        suptitle('Train on CC PYR');
                    elseif i==2
                        suptitle('Train on Alt1 PYR');
                    end
                end
                
                % other muscles
                figure;
                subplot(2,3,1); hold on;
                plot(Y_real_input{1,1}(:,3),'k');
                ylabel('pal');
                title('cc');
                plot(Y_fit_input{1,1}(:,3),'g');
                axis square;
                box off;
                
                subplot(2,3,2); hold on;
                plot(Y_real_input{1,2}(:,3),'k');
                title('alt1');
                plot(Y_fit_input{1,2}(:,3),'g');
                axis square;
                box off;
                
                subplot(2,3,3); hold on;
                plot(Y_real_input{1,3}(:,3),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,3),'g');
                axis square;
                box off;

                subplot(2,3,4); hold on;
                plot(Y_real_input{1,1}(:,2),'k');
                ylabel('edc');
                plot(Y_fit_input{1,1}(:,2),'g');
                axis square;
                box off;

                subplot(2,3,5); hold on;
                plot(Y_real_input{1,2}(:,2),'k');
                plot(Y_fit_input{1,2}(:,2),'g');
                axis square;
                box off;
                
                subplot(2,3,6); hold on;
                plot(Y_real_input{1,3}(:,2),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,2),'g');
                axis square;
                box off;
                
                if k == 1
                    if i==1
                        suptitle('Train on CC ALL');
                    elseif i==2
                        suptitle('Train on Alt1 ALL');
                    end
                elseif k == 2
                    if i==1
                        suptitle('Train on CC PYR');
                    elseif i==2
                        suptitle('Train on Alt1 PYR');
                    end
                end

                figure;
                subplot(2,3,1); hold on;
                plot(Y_real_input{1,1}(:,4),'k');
                ylabel('pec');
                title('cc');
                plot(Y_fit_input{1,1}(:,4),'g');
                axis square;
                box off;
                
                subplot(2,3,2); hold on;
                plot(Y_real_input{1,2}(:,4),'k');
                title('alt1');
                plot(Y_fit_input{1,2}(:,4),'g');
                axis square;
                box off;
                
                subplot(2,3,3); hold on;
                plot(Y_real_input{1,3}(:,4),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,4),'g');
                axis square;
                box off;
                
                subplot(2,3,4); hold on;
                plot(Y_real_input{1,1}(:,5),'k');
                ylabel('spdel');
                plot(Y_fit_input{1,1}(:,5),'g');
                axis square;
                box off;
                
                subplot(2,3,5); hold on;
                plot(Y_real_input{1,2}(:,5),'k');
                plot(Y_fit_input{1,2}(:,5),'g');
                axis square;
                box off;
                
                subplot(2,3,6); hold on;
                plot(Y_real_input{1,3}(:,5),'k');
                title('alt2');
                plot(Y_fit_input{1,3}(:,5),'g');
                axis square;
                box off;
                
                if k == 1
                    if i==1
                        suptitle('Train on CC ALL');
                    elseif i==2
                        suptitle('Train on Alt1 ALL');
                    end
                elseif k == 2
                    if i==1
                        suptitle('Train on CC PYR');
                    elseif i==2
                        suptitle('Train on Alt1 PYR');
                        
                    end
                end
            end
            
        end
        
        
        if k == 1
            Y_real_CCmodel_all=Y_real_CCmodel;
            Y_fit_CCmodel_all=Y_fit_CCmodel;
            M_CCmodel_all=M_CCmodel;
            Test_R2_CCmodel_all=Test_R2_CCmodel;
            lambda_opt_CCmodel_all=lambda_opt_CCmodel;
            Y_fit_CCmodel_Xperms_all=Y_fit_CCmodel_Xperms;
            
            Y_real_A1model_all=Y_real_A1model;
            Y_fit_A1model_all=Y_fit_A1model;
            M_A1model_all=M_A1model;
            Test_R2_A1model_all=Test_R2_A1model;
            lambda_opt_A1model_all=lambda_opt_A1model;
            Y_fit_A1model_Xperms_all=Y_fit_A1model_Xperms;
            
            M_A2model_all=M_A2model;
            
        elseif k == 2
            Y_real_CCmodel_pyr=Y_real_CCmodel;
            Y_fit_CCmodel_pyr=Y_fit_CCmodel;
            M_CCmodel_pyr=M_CCmodel;
            Test_R2_CCmodel_pyr=Test_R2_CCmodel;
            lambda_opt_CCmodel_pyr=lambda_opt_CCmodel;
            Y_fit_CCmodel_Xperms_pyr=Y_fit_CCmodel_Xperms;
            
            Y_real_A1model_pyr=Y_real_A1model;
            Y_fit_A1model_pyr=Y_fit_A1model;
            M_A1model_pyr=M_A1model;
            Test_R2_A1model_pyr=Test_R2_A1model;
            lambda_opt_A1model_pyr=lambda_opt_A1model;
            Y_fit_A1model_Xperms_pyr=Y_fit_A1model_Xperms;
            
            M_A2model_pyr=M_A2model;
        end
    end
    
    
    
    %% linear regression of joint model naive lambda
    warning('off','stats:regress:RankDefDesignMat');
    
    for k = 1:2 %1 for all neurons, 2 for just pyr
        if k == 1
            input=Data;
            disp('lin reg joint all');
        elseif k == 2
            input=Data_pyr;
            disp('lin reg joint pyr');
        end
        R2_conds=nan*ones(length(conds));
        M_opt=cell(length(conds),1);
        Y_train_opt=cell(length(conds),1);
        Y=cell(length(conds),length(conds));
        Y_fit=cell(length(conds),length(conds));
        
        Train_R2=nan*ones(length(lambda_vec),1);
        Val_R2=nan*ones(length(lambda_vec),1);
        M_local=cell(length(lambda_vec));
        Test_R2=nan*ones(num_perms_lambda,length(conds));
        Y_fit=cell(num_perms_lambda,length(conds));
        Y_real=cell(num_perms_lambda,length(conds));
        lambda_opt=nan*ones(1,num_perms_lambda);
        M=cell(num_perms_lambda);
        
        % here, lambda calculated from train and val and tested on test trials in EACH permutation
        for b=1:num_perms_lambda
            %b
            
            % Fit the model on the training set of Alt1-Cc concatenated in time
            % neural data from Alt1 and Cc concat in time, mean center together
            X_train=[input(2).train_neur{b}; input(1).train_neur{b}];
            X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
            
            %muscle data from Alt1 and Cc concat in time, mean center together
            Y_train=[input(2).train_emg{b}; input(1).train_emg{b}];
            Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
            
            % neural data from Alt1 and Cc concat in time, mean center together
            X_val=[input(2).val_neur{b}; input(1).val_neur{b}];
            X_val=X_val-repmat(mean(X_val,1),size(X_val,1),1); % mean center X_val
            
            % muscle data from Alt1 and Cc concat in time, mean center together
            Y_val=[input(2).val_emg{b}; input(1).val_emg{b}];
            Y_val=Y_val-repmat(mean(Y_val,1),size(Y_val,1),1); % mean center Y_val

            for j=1:length(lambda_vec)
                [Train_R2(j), M_local{j}, ~] = anyRegress( X_train, Y_train, 'ridge',10^lambda_vec(j));
                [Val_R2(j),~,~] = anyRegress( X_val, Y_val, 'M',M_local{j} );
            end
            % pick best lambda
            [max_R2,lambda_opt_ind]=max(Val_R2); %find lambda that gave highest R2
            lambda_opt(b)=10^(lambda_vec(lambda_opt_ind));
            
            %now combine the test and val trials
            X_train=([input(2).train_neur{b}; input(1).train_neur{b}]+[input(2).val_neur{b}; input(1).val_neur{b}])./2; %mean of test and val
            X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
            Y_train=([input(2).train_emg{b}; input(1).train_emg{b}]+[input(2).val_emg{b}; input(1).val_emg{b}])./2;%mean of test and val
            Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
            
            
            % get the M for train data, use within-iteration lambda
            [~, M{b}, ~] = anyRegress(X_train, Y_train, 'ridge',lambda_opt(b));
            
            % test the weights
            for c=1:length(conds)
                % use only subset of trials (test)
                X=input(c).test_neur{b};
                X=X-repmat(mean(X,1),size(X,1),1); % mean center X_val
                Y=[input(c).test_emg{b}];
                Y=Y-repmat(mean(Y,1),size(Y,1),1); % mean center Y_train
                Y_real{b,c}=Y; %record your real emg for examples
                [Test_R2(b,c),~,Y_fit{b,c}]=anyRegress(X, Y, 'M',M{b} );
            end
            
            
        end
        
        %average over permutations
        Test_R2_mean=mean(Test_R2); %average over permutations
        %identify the iteration whose R2 vals were closest in aggregate to mean
        Test_R2_diff=nan*ones(size(Test_R2));
        for j=1:size(Test_R2,1)
            Test_R2_diff(j,:)=abs(Test_R2(j,:)-Test_R2_mean);
        end
        Test_R2_diff_sum=sum(Test_R2_diff(:,1:3),2);
        [val,ind]=min(Test_R2_diff_sum);
        
        Y_real_joint=Y_real(ind,:);
        Y_fit_joint=Y_fit(ind,:);
        M_joint=M;
        Test_R2_joint=Test_R2_mean;
        lambda_opt_joint=lambda_opt;
        
        figure;
        bar(Test_R2_joint(1:3));
        set(gca,'fontsize',20)
        xlabel('Condition used for fitting M');
        xticklabels(conds);
        ax=gca;
        set(ax,'XTickLabel',conds);
        legend(conds(1:3));
        %ylim([-2 1]);
        if k == 1
            title('R^2 ALL neurons');
        elseif k == 2
            title('R^2 PYR neurons');
        end
        
        if plot_lin_reg_traces
            ylim_vals=[-0.2 0.26];
            Y_real_input=Y_real_joint;
            Y_fit_input=Y_fit_joint;
            
            figure;
            subplot(2,3,1); hold on;
            plot(Y_real_input{1,1}(:,1),'k');
            ylabel('bi');
            title('cc');
            plot(Y_fit_input{1,1}(:,1),'g');
            axis square;
            box off;
            
            subplot(2,3,2); hold on;
            plot(Y_real_input{1,2}(:,1),'k');
            title('alt1');
            plot(Y_fit_input{1,2}(:,1),'g');
            axis square;
            box off;
            
            subplot(2,3,3); hold on;
            plot(Y_real_input{1,3}(:,1),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,1),'g');
            axis square;
            box off;

            subplot(2,3,4); hold on;
            plot(Y_real_input{1,1}(:,6),'k');
            ylabel('tri');
            plot(Y_fit_input{1,1}(:,6),'g');
            axis square;
            box off;
            
            subplot(2,3,5); hold on;
            plot(Y_real_input{1,2}(:,6),'k');
            plot(Y_fit_input{1,2}(:,6),'g');
            axis square;
            box off;
            
            subplot(2,3,6); hold on;
            plot(Y_real_input{1,3}(:,6),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,6),'g');
            axis square;
            box off;
            
            if k == 1
                if i==1
                    suptitle('Train on CC ALL');
                elseif i==2
                    suptitle('Train on Alt1 ALL');
                end
            elseif k == 2
                if i==1
                    suptitle('Train on CC PYR');
                elseif i==2
                    suptitle('Train on Alt1 PYR');
                end
            end
            
            % other muscles
            figure;
            subplot(2,3,1); hold on;
            plot(Y_real_input{1,1}(:,3),'k');
            ylabel('pal');
            title('cc');
            plot(Y_fit_input{1,1}(:,3),'g');
            axis square;
            box off;
            
            subplot(2,3,2); hold on;
            plot(Y_real_input{1,2}(:,3),'k');
            title('alt1');
            plot(Y_fit_input{1,2}(:,3),'g');
            axis square;
            box off;
            
            subplot(2,3,3); hold on;
            plot(Y_real_input{1,3}(:,3),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,3),'g');
            axis square;
            box off;
            
            subplot(2,3,4); hold on;
            plot(Y_real_input{1,1}(:,2),'k');
            ylabel('edc');
            plot(Y_fit_input{1,1}(:,2),'g');
            axis square;
            box off;
            
            subplot(2,3,5); hold on;
            plot(Y_real_input{1,2}(:,2),'k');
            plot(Y_fit_input{1,2}(:,2),'g');
            axis square;
            box off;
            
            subplot(2,3,6); hold on;
            plot(Y_real_input{1,3}(:,2),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,2),'g');
            axis square;
            box off;
            
            if k == 1
                if i==1
                    suptitle('Train on CC ALL');
                elseif i==2
                    suptitle('Train on Alt1 ALL');
                end
            elseif k == 2
                if i==1
                    suptitle('Train on CC PYR');
                elseif i==2
                    suptitle('Train on Alt1 PYR');
                end
            end
            
            figure;
            subplot(2,3,1); hold on;
            plot(Y_real_input{1,1}(:,4),'k');
            ylabel('pec');
            title('cc');
            plot(Y_fit_input{1,1}(:,4),'g');
            axis square;
            box off;

            subplot(2,3,2); hold on;
            plot(Y_real_input{1,2}(:,4),'k');
            title('alt1');
            plot(Y_fit_input{1,2}(:,4),'g');
            axis square;
            box off;
            
            subplot(2,3,3); hold on;
            plot(Y_real_input{1,3}(:,4),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,4),'g');
            axis square;
            box off;
            
            subplot(2,3,4); hold on;
            plot(Y_real_input{1,1}(:,5),'k');
            ylabel('spdel');
            plot(Y_fit_input{1,1}(:,5),'g');
            axis square;
            box off;
            
            subplot(2,3,5); hold on;
            plot(Y_real_input{1,2}(:,5),'k');
            plot(Y_fit_input{1,2}(:,5),'g');
            axis square;
            box off;
            
            subplot(2,3,6); hold on;
            plot(Y_real_input{1,3}(:,5),'k');
            title('alt2');
            plot(Y_fit_input{1,3}(:,5),'g');
            axis square;
            box off;
            
            if k == 1
                suptitle('Train on joint ALL');
            elseif k == 2
                suptitle('Train on joint PYR');
                
            end
        end

        
        if k == 1
            Y_real_joint_all=Y_real_joint;
            Y_fit_joint_all=Y_fit_joint;
            M_joint_all=M_joint;
            Test_R2_joint_all=Test_R2_joint;
            lambda_opt_joint_all=lambda_opt_joint;
        elseif k == 2
            Y_real_joint_pyr=Y_real_joint;
            Y_fit_joint_pyr=Y_fit_joint;
            M_joint_pyr=M_joint;
            Test_R2_joint_pyr=Test_R2_joint;
            lambda_opt_joint_pyr=lambda_opt_joint;
        end
    end
    
    
    %% joint model control: scale the cocon vals on neural and emg
    
    scale_factors=[0.1:0.1:1.3];
    warning('off','stats:regress:RankDefDesignMat');
    
    for k = 1:2 %1 for all neurons, 2 for just pyr
        if k == 1
            input=Data;
            disp('lin reg joint scale all');
        elseif k == 2
            input=Data_pyr;
            disp('lin reg joint scale pyr');
        end
        R2_conds_scale=nan*ones(length(scale_factors),length(conds));
        for scale_it=1:length(scale_factors)
            R2_conds=nan*ones(length(conds));
            M_opt=cell(length(conds),1);
            Y_train_opt=cell(length(conds),1);
            Y=cell(length(conds),length(conds));
            Y_fit=cell(length(conds),length(conds));
            
            Train_R2=nan*ones(length(lambda_vec),1);
            Val_R2=nan*ones(length(lambda_vec),1);
            M_local=cell(length(lambda_vec));
            % new
            Test_R2=nan*ones(num_perms_lambda,length(conds));
            Y_fit=cell(num_perms_lambda,length(conds));
            Y_real=cell(num_perms_lambda,length(conds));
            lambda_opt=nan*ones(1,num_perms_lambda);
            M=cell(num_perms_lambda);
            
            % here, lambda calculated from train and val and tested on test trials in EACH permutation
            for b=1:num_perms_lambda
                %b
                % Fit the model on the training set of Alt1-Cc concatenated in time
                % neural data from Alt1 and Cc concat in time, mean center together
                X_train=[input(2).train_neur{b}; input(1).train_neur{b}*scale_factors(scale_it)];
                X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
                
                %muscle data from Alt1 and Cc concat in time, mean center together
                Y_train=[input(2).train_emg{b}; input(1).train_emg{b}*scale_factors(scale_it)];
                Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
                
                % neural data from Alt1 and Cc concat in time, mean center together
                X_val=[input(2).val_neur{b}; input(1).val_neur{b}*scale_factors(scale_it)];
                X_val=X_val-repmat(mean(X_val,1),size(X_val,1),1); % mean center X_val
                
                % muscle data from Alt1 and Cc concat in time, mean center together
                Y_val=[input(2).val_emg{b}; input(1).val_emg{b}*scale_factors(scale_it)];
                Y_val=Y_val-repmat(mean(Y_val,1),size(Y_val,1),1); % mean center Y_val
                
                for j=1:length(lambda_vec)
                    [Train_R2(j), M_local{j}, ~] = anyRegress( X_train, Y_train, 'ridge',10^lambda_vec(j));
                    [Val_R2(j),~,~] = anyRegress( X_val, Y_val, 'M',M_local{j} );
                end
                % pick best lambda
                [max_R2,lambda_opt_ind]=max(Val_R2); %find lambda that gave highest R2
                lambda_opt(b)=10^(lambda_vec(lambda_opt_ind)); %changed bc we're using powers of 10
                
                %now combine the test and val trials
                X_train=([input(2).train_neur{b}; input(1).train_neur{b}*scale_factors(scale_it)]+[input(2).val_neur{b}; input(1).val_neur{b}*scale_factors(scale_it)])./2; %mean of test and val
                X_train=X_train-repmat(mean(X_train,1),size(X_train,1),1); % mean center X_train
                Y_train=([input(2).train_emg{b}; input(1).train_emg{b}*scale_factors(scale_it)]+[input(2).val_emg{b}; input(1).val_emg{b}*scale_factors(scale_it)])./2;%mean of test and val
                Y_train=Y_train-repmat(mean(Y_train,1),size(Y_train,1),1); % mean center Y_train
                
                % get the M for train data, use within-iteration lambda
                [~, M{b}, ~] = anyRegress(X_train, Y_train, 'ridge',lambda_opt(b));
                
                % test the weights
                for c=1:length(conds)
                    if c==1
                        X=input(c).test_neur{b}*scale_factors(scale_it);
                        Y=input(c).test_emg{b}*scale_factors(scale_it);
                    else
                        X=input(c).test_neur{b};
                        Y=[input(c).test_emg{b}];
                    end

                    X=X-repmat(mean(X,1),size(X,1),1); % mean center X_val
                    Y=Y-repmat(mean(Y,1),size(Y,1),1); % mean center Y_train
                    Y_real{b,c}=Y; %record your real emg for examples
                    [Test_R2(b,c),~,Y_fit{b,c}]=anyRegress(X, Y, 'M',M{b} );
                end
            end
            
            %average over permutations
            R2_conds_scale(scale_it,:)=mean(Test_R2);
        end
        
        figure;
        hold on;
        plot(scale_factors,R2_conds_scale(:,1),'r');
        plot(scale_factors,R2_conds_scale(:,2),'b');
        plot(scale_factors,R2_conds_scale(:,3),'k');
        legend('cc','a1','a2');
        xlabel('fraction of CC signal');
        ylabel('R^2');
        
        if k == 1
            title('Scaled CC signals ALL');
            R2_conds_joint_scale_all=R2_conds_scale;
        elseif k == 2
            title('Scaled CC signals PYR');
            R2_conds_joint_scale_pyr=R2_conds_scale;
        end
        
    end
    
    
    
    %% Extract muscle-related neural signal using model

    for k = 1:2 %1 for all neurons, 2 for just pyr
        if k == 1
            input=Data;
            disp('cleaning all');
        elseif k == 2
            input=Data_pyr;
            disp('cleaning pyr');
        end
        VarCap=zeros(6,4);
        for ii=1:length(conds)
            if ii==4
                if k == 1
                    Data(ii).clean_data=input(ii).fr;
                elseif k == 2
                    Data_pyr(ii).clean_data=input(ii).fr;
                end
            else
            
            %1) use the mean of previous permutations of weights for these weights
            if k == 1
                if ii == 1
                    w_perms=M_CCmodel_all;
                elseif ii == 2
                    w_perms=M_A1model_all;
                elseif ii == 3
                    w_perms=M_A2model_all;
                end
            elseif k == 2
                if ii == 1
                    w_perms=M_CCmodel_pyr;
                elseif ii == 2
                    w_perms=M_A1model_pyr;
                elseif ii == 3
                    w_perms=M_A2model_pyr;
                end
            end
            w_perms(:,2)=[];
            [row,col]=size(w_perms);
            mean_w=w_perms{1}; %initialize with 1st perm's weights
            for iter=2:row
                mean_w=mean_w+w_perms{iter};
            end
            weights=mean_w/row;

            %2) orthonormailze these weight vectors using graham-schmidt process
            orthoWeights = Gram_Schmidt_Process(weights);
            
            %3) determine variance captured
            neural_data=input(ii).fr-repmat(mean(input(ii).fr),size(input(ii).fr,1),1); %mean-centered
            for i=1:6
                fitsData = neural_data*orthoWeights(:,i);
                VarCap(i,ii) = trace(cov(fitsData)) / trace(cov(neural_data));
            end
            
            %3b determine the n weights that carry most var
            [valsSort,indSort]=sort(VarCap(:,ii));
            weights_select=indSort(end-(numPCs_forClean-1):end);
            
            %4) Project TxN data onto Nxn matrix of n orthonormal vectors
            proj_data=neural_data*orthoWeights(:,weights_select);
            
            %5)Project out of orthonormal subspace with transpose of orthonormal vectors
            clean_data=proj_data*orthoWeights(:,weights_select)';
            
            %add back mean that had been subtracted
            clean_data=clean_data+repmat(mean(input(ii).fr),size(input(ii).fr,1),1);
            neural_data=neural_data+repmat(mean(input(ii).fr),size(input(ii).fr,1),1);
            
            %                 if plot_cleaned_neuron_examples
            %                     %6) plot some example cell timeseries pre- and post-cleaning to compare
            %                     cellnum=[1,5,10,15,20,25,30,35];
            %                     figure;
            %                     for i=1:8
            %                         subplot(4,2,i)
            %                         hold on;
            %                         plot(neural_data(:,cellnum(i)),'k');
            %                         plot(clean_data(:,cellnum(i)),'g');
            %                     end
            %                     legend('Original','Cleaned');
            %                     if k == 1
            %                         suptitle(strcat('All neurons-',conds{ii}));
            %                     elseif k == 2
            %                         suptitle(strcat('Pyr neurons-',conds{ii}));
            %                     end
            %                 end
            if k == 1
                Data(ii).clean_data=clean_data;
            elseif k == 2
                Data_pyr(ii).clean_data=clean_data;
            end
            end
            
        end
    end
    
    
    
    %% PCA
    
    for kk=1:2
        if kk==1
            input=Data;
        elseif kk==2
            input=Data_pyr;
        end
        warning('off','stats:pca:ColRankDefX')
        totalvar=zeros(length(conds),1);
        average_dist=zeros(length(conds),1);
        projs=cell(length(conds),1);
        PCs=cell(length(conds),1);
        explained_mat=zeros(length(conds),numPCs);
        
        % new
        all_mat=vertcat(input(:).fr);
        mean_across_conds=mean(all_mat);
        
        for i=1:length(conds)
            fr_mat=input(i).fr-repmat(mean_across_conds,size(input(i).fr,1),1);

            [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
            PCs{i}=PC;
            projs{i}=scores(:,1:numPCs);
            eigs=sort(eig(cov(scores)),'descend');
            totalvar(i)=sum(eig(cov(scores)));
            average_dist(i)=(mean((max(abs(fr_mat),[],1))));
            explained_mat(i,:)=explaineds(1:numPCs);
        end
        figure;
        co=colormap(jet(length(conds)));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        plot(cumsum(explained_mat'),'linewidth',lw);
        legend(conds{1:4});
        if kk==1
            title('Variance Explained (Neural) ALL NEURONS');
        elseif kk==2
            title('Variance Explained (Neural) PYR');
        end
        set(gca,'fontsize',14);
        ylim([min(min(cumsum(explained_mat')))*0.95 100]);
        if ~isempty(totalvar)
            axes('Position',[.45 .25 .3 .3])
            box on
            plot(1:length(conds),totalvar,'-k','linewidth',lw);
            ax=gca;
            set(ax,'XTick',1:length(conds));
            set(ax,'XTickLabel',{'CC','Alt1','Alt2','Inact'});
            xlabel('Condition'); ylabel('Total Variance');set(gca,'fontsize',14);
            ylim([min(totalvar)*0.95 max(totalvar)*1.05]);
        end
        
        % Analyze variance explained by PCs of different conditions, INCLUDING INACT WITH PC breakdown
        explaineds=cell(length(conds),1);
        condfig=nan(length(conds),1);
        for i=1:length(conds)
            explaineds{i}=zeros(length(conds),numPCs);
            % Project each forward condition into alt and cc PCs
            condfig(i)=figure;
            co=colormap(lines(length(conds)));
            set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
            PCs_cond=PCs{i};
            for c=1:length(conds)
                fr_mat=input(c).fr-repmat(mean_across_conds,size(input(c).fr,1),1);
                
                for cc=1:numPCs
                    scores_thiscond=fr_mat*PCs_cond(:,1:cc);
                    eigs_thiscond=eig(cov(scores_thiscond));
                    explaineds{i}(c,cc)=sum(eigs_thiscond)/totalvar(c);
                end
                figure(condfig(i));
                if numPCs>2
                    plot3(scores_thiscond(:,1),scores_thiscond(:,2),scores_thiscond(:,3),'linewidth',lw); hold on;
                elseif numPCs>1
                    plot(scores_thiscond(:,1),scores_thiscond(:,2),'linewidth',lw); hold on;
                end
                fprintf('Percent variance explained by %s for %s : %d \n',conds{i},conds{c},round(100*explaineds{i}(c)))
            end
            figure(condfig(i));
            set(gca,'fontsize',18);
            
            if kk==1
                title(strcat('All Cons proj in ',{' '},conds{i},{' '},' PC space ALL NEURONS'));
            elseif kk==2
                title(strcat('All Cons proj in ',{' '},conds{i},{' '},' PC space PYR'));
            end
            
            xlabel(strcat(conds{i},{' '},' PC1')); ylabel(strcat(conds{i},{' '},' PC2')); zlabel(strcat(conds{i},{' '},' PC3'));
            legend(conds);
        end
        
        if kk==1
            explaineds_all=explaineds;
        elseif kk==2
            explaineds_pyr=explaineds;
        end
        
        %make plottable ver of cell array
        explaineds_2=explaineds;
        for i=1:length(explaineds)
            explaineds_2{i}=explaineds{i}(:,end);
        end
        
        %         figure;
        %         bar(1:length(conds),horzcat(explaineds_2{:})');
        %         set(gca,'fontsize',18);
        %         legend(strcat(conds,' PCs'));
        %         ylim([0 1]);
        %         xlabel('Conditions');
        %         %xticklabels(conds);
        %         ax=gca; set(ax,'XTickLabel',conds);
        %         if kk==1
        %             title('Var explained in Diff Conds ALL NEURONS');
        %         elseif kk==2
        %             title('Var explained in Diff Conds PYR');
        %     %     elseif kk==3
        %     %         title('Var explained in Diff Conds INs');
        %         end
        
        
        PCsVarExplained=[{'CC';'Alt1';'Alt2'} explaineds_2(1:3)];
        PCsVarExplained_vals=PCsVarExplained(:,2);
        
        figure;
        bar([1 2 4 5 7 8],[PCsVarExplained_vals{1}(2)/PCsVarExplained_vals{1}(1) PCsVarExplained_vals{1}(3)/PCsVarExplained_vals{1}(1)...
            PCsVarExplained_vals{2}(1)/PCsVarExplained_vals{2}(2) PCsVarExplained_vals{2}(3)/PCsVarExplained_vals{2}(2)...
            PCsVarExplained_vals{3}(1)/PCsVarExplained_vals{3}(3) PCsVarExplained_vals{3}(2)/PCsVarExplained_vals{3}(3)]);
        ax=gca; set(ax,'XTickLabel',{'a1/cc','a2/cc','cc/a1','a2/a1','cc/a2','a1/a2'});
        ylim([0 1]);
        text(0.6,0.9,num2str([PCsVarExplained_vals{1}(2)/PCsVarExplained_vals{1}(1) PCsVarExplained_vals{1}(3)/PCsVarExplained_vals{1}(1)...
            PCsVarExplained_vals{2}(1)/PCsVarExplained_vals{2}(2) PCsVarExplained_vals{2}(3)/PCsVarExplained_vals{2}(2)...
            PCsVarExplained_vals{3}(1)/PCsVarExplained_vals{3}(3) PCsVarExplained_vals{3}(2)/PCsVarExplained_vals{3}(3)]));
        
        varVals=[PCsVarExplained_vals{1}(2)/PCsVarExplained_vals{1}(1) PCsVarExplained_vals{1}(3)/PCsVarExplained_vals{1}(1)...
            PCsVarExplained_vals{2}(1)/PCsVarExplained_vals{2}(2) PCsVarExplained_vals{2}(3)/PCsVarExplained_vals{2}(2)...
            PCsVarExplained_vals{3}(1)/PCsVarExplained_vals{3}(3) PCsVarExplained_vals{3}(2)/PCsVarExplained_vals{3}(3)];
        
        if kk==1
            title('ALL NEURONS');
        elseif kk==2
            title('PYR');
        end
        
        %var over PCs (cumsum)
        figure;
        subplot(1,4,1);
        plot(explaineds{1}');
        ylim([0 1]);
        title('CC PCs');
        subplot(1,4,2);
        plot(explaineds{2}');
        ylim([0 1]);
        title('A1 PCs');
        subplot(1,4,3);
        plot(explaineds{3}');
        ylim([0 1]);
        title('A2 PCs');
        subplot(1,4,4);
        plot(explaineds{4}');
        ylim([0 1]);
        legend('cc','a1','a2','inact');
        title('Inact PCs');
        if kk==1
            suptitle('ALL NEURONS');
        elseif kk==2
            suptitle('PYRs');
        end
        
        
        
        % PCA for all directions together, INCLUDING INACT
        % find the mean across all 4 conds, and subtract that off, but only project in the 2 conds (A1 and CC)
        
        allA=vertcat(input.fr); %get mean with 4 conds to center
        meanA=mean(allA,1);
        allA=vertcat(input([1:2]).fr); %make PC space w/o Alt2 or Inact
        allA = allA-repmat(meanA,size(allA,1),1);
        [PCs,scores] = pca(allA);
        figure;
        co=colormap(lines(length(conds)));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        for i=1:length(conds)
            if size(meanA)==size(input(i).fr)
                fr_mat=input(i).fr-meanA;
            else
                fr_mat=input(i).fr-repmat(meanA,size(input(i).fr,1),1);
            end
            scores=fr_mat*PCs(:,1:numPCs);
            if numPCs>2
                plot3(scores(:,1),scores(:,2),scores(:,3),'linewidth',lw); hold on;
            elseif numPCs>1
                plot(scores(:,1),scores(:,2),'linewidth',lw); hold on;
            end
        end
        set(gca,'fontsize',18);
        
        if kk==1
            title('All Conds proj in Alt1-CC PC space ALL NEURONS');
            if mouse==1
                xlim([-60 65]);
                ylim([-60 60]);
                zlim([-60 60]);
                xticks([-60 0 60])
                yticks([-60 0 60])
                zticks([-60 0 60])
            end
        elseif kk==2
            title('All Conds proj in Alt1-CC PC space PYR');
            if mouse==1
                xlim([-60 60.1]);
                ylim([-60 60]);
                zlim([-60 60]);
                xticks([-60 0 60])
                yticks([-60 0 60])
                zticks([-60 0 60])
            end
        end
        legend(conds);
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        grid on;
        %rotate3d on;
    end

 
    %% ratio of variance

    % find top 6 PCs for [Alt1; Alt2; Cocon; Cocon]
    % do PCA, keep PCs, throw out "scores"
    % PCs will be in neuron num x neuron num matrix
    % project the data in these PCs (matrix V): (the PCs common to all epochs)
    % matrix A (conds x 300) * V(:,1:6) = A' ct X 6
    
    % find top 6 PCs for [Alt1; Alt2; Cocon; Cocon]
    input=Data;
    warning('off','stats:pca:ColRankDefX')
    totalvar=zeros(1,1);
    average_dist=zeros(1,1);
    projs=cell(1,1);
    PCs=cell(1,1);
    explained_mat=zeros(1,numPCs);
    
    % make matrix of Alt1; Alt2; Cc; Cc
    fr_mat=[input(2).fr; input(3).fr; input(1).fr; input(1).fr];
    % mean center, only subtract off total mean for each neuron across the conditions
    fr_mat=fr_mat-repmat(mean(fr_mat),size(fr_mat,1),1);
    
    [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
    PCs=PC;
    projs{i}=scores(:,1:numPCs);
    eigs=sort(eig(cov(scores)),'descend');
    totalvar(i)=sum(eig(cov(scores)));
    average_dist(i)=(mean((max(abs(fr_mat),[],1))));
    explained_mat(i,:)=explaineds(1:numPCs);
    
    % project the data in these PCs (matrix V):
    % matrix A (conds x 300) * V(:,1:6) = A' ct X 6
    % take off 2nd CC from fr_mat
    fr_mat=fr_mat(1:900,:);
    
    A_prime = fr_mat*PCs(:,1:6);
    % mean center across all conditions
    A_prime_centered=A_prime-repmat(mean(A_prime),size(A_prime,1),1);
    
    a1_condense=A_prime_centered(1:300,:);
    a2_condense=A_prime_centered(301:600,:);
    cc_condense=A_prime_centered(601:900,:);
    
    % now do PCA on these data
    figure;
    subplot(1,2,1);
    hold on;
    
    % neural
    warning('off','stats:pca:ColRankDefX')
    totalvar=zeros(3,1);
    average_dist=zeros(3,1);
    projs=cell(3,1);
    PCs=cell(3,1);
    explained_mat=zeros(3,numPCs);
    for i=1:3
        if i==1
            fr_mat=cc_condense;
        elseif i==2
            fr_mat=a1_condense;
        elseif i==3
            fr_mat=a2_condense;
        end
        [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
        PCs{i}=PC;
        projs{i}=scores(:,1:numPCs);
        eigs=sort(eig(cov(scores)),'descend');
        totalvar(i)=sum(eig(cov(scores)));
        average_dist(i)=(mean((max(abs(fr_mat),[],1))));
        explained_mat(i,:)=explaineds(1:numPCs);
    end
    
    % epoch by epoch
    plot(cumsum(explained_mat(1,:)'),'r','linewidth',lw);
    plot(cumsum(explained_mat(2,:)'),'b','linewidth',lw);
    plot(cumsum(explained_mat(3,:)'),'g','linewidth',lw);
    set(gca,'fontsize',14);
    ylim([50 100]);
    hold on;
    
    % mean across epochs
    subplot(1,2,2);
    hold on;
    plot(mean([cumsum(explained_mat(1,:)') cumsum(explained_mat(2,:)') cumsum(explained_mat(3,:)')]'),'k','linewidth',lw);
    set(gca,'fontsize',14);
    ylim([50 100]);
    hold on;
    
    ortho_control_neuro=mean([cumsum(explained_mat(1,:)') cumsum(explained_mat(2,:)') cumsum(explained_mat(3,:)')]');
    
    
    % add emg in dotted line
    warning('off','stats:pca:ColRankDefX')
    totalvar=zeros(length(conds),1);
    average_dist=zeros(length(conds),1);
    projs=cell(length(conds),1);
    PCs=cell(length(conds),1);
    explained_mat=zeros(length(conds),numPCs_emg);
    
    % make matrix of muscles across conditions
    x=vertcat(Data(1:3).emg);
    
    mean_EMG_across_conds=mean(vertcat(Data(1:3).emg));
    
    for i=1:length(conds)-1
        fr_mat=Data(i).emg-repmat(mean_EMG_across_conds,size(Data(i).emg,1),1);
        [PC,scores,latents,coeffs,explaineds_emg]=pca(fr_mat);
        PCs{i}=PC;
        projs{i}=scores(:,1:numPCs_emg);
        eigs=sort(eig(cov(scores)),'descend');
        totalvar(i)=sum(eig(cov(scores)));
        average_dist(i)=(mean((max(abs(fr_mat),[],1))));
        explained_mat(i,:)=explaineds_emg(1:numPCs_emg);
    end
    
    subplot(1,2,1);
    hold on;
    plot(cumsum(explained_mat(1,:)'),'r','linewidth',lw,'linestyle','--');
    plot(cumsum(explained_mat(2,:)'),'b','linewidth',lw,'linestyle','--');
    plot(cumsum(explained_mat(3,:)'),'g','linewidth',lw,'linestyle','--');
    legend({'neural cc','neural alt1','neural alt2','emg cc','emg alt1','emg alt2'});
    set(gca,'fontsize',14);
    ylim([55 101]);
    xlim([0.8 6.2]);
    ylabel('Cumulative variance explained');
    xlabel('PCs');
    axis square;
    box off;
    
    subplot(1,2,2);
    hold on;
    plot(mean([cumsum(explained_mat(1,:)') cumsum(explained_mat(2,:)') cumsum(explained_mat(3,:)')]'),'k','linewidth',lw,'linestyle','--');
    legend({'neural mean','emg mean'});
    set(gca,'fontsize',14);
    ylim([55 101]);
    xlim([0.8 6.2]);
    ylabel('Cumulative variance explained');
    xlabel('PCs');
    axis square;
    box off;
    suptitle('Neural data compressed to 6 PCs and EMG');
    
    ortho_control_emg=mean([cumsum(explained_mat(1,:)') cumsum(explained_mat(2,:)') cumsum(explained_mat(3,:)')]');
    
    
    % now you want to use the "conditioned/compressed" data but find the PCs
    % that are epoch specific and test var ex for each epoch in these
    
    % ultimately you want:
    % Variance of co-contraction data captured by top two (or three) alternation PCs.
    % Variance of co-contraction data captured by top two (or three) co-contraction PCs.
    % Ratio of the above. Compare that ratio for neural and EMG data.
    
    ylim_max=max(max([a1_condense; a2_condense; cc_condense]));
    ylim_min=min(min([a1_condense; a2_condense; cc_condense]));
    
    figure;
    subplot(2,3,1);
    plot(a1_condense(:,1:3));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    subplot(2,3,2);
    plot(a2_condense(:,1:3));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    subplot(2,3,3);
    plot(cc_condense(:,1:3));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    legend({'PC1','PC2','PC3'});
    
    subplot(2,3,4);
    plot(a1_condense(:,4:6));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    subplot(2,3,5);
    plot(a2_condense(:,4:6));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    subplot(2,3,6);
    plot(cc_condense(:,4:6));
    ylim([ylim_min ylim_max]);
    axis square;
    box off;
    legend({'PC4','PC5','PC6'});
    
    num_conds=3; % no inact data currently
    
    warning('off','stats:pca:ColRankDefX')
    totalvar=zeros(num_conds,1);
    average_dist=zeros(num_conds,1);
    projs=cell(num_conds,1);
    PCs=cell(num_conds,1);
    explained_mat=zeros(num_conds,numPCs);
    for i=1:num_conds
        if i==1
            fr_mat=cc_condense;
        elseif i==2
            fr_mat=a1_condense;
        elseif i==3
            fr_mat=a2_condense;
        end
        
        [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
        PCs{i}=PC;
        projs{i}=scores(:,1:numPCs);
        eigs=sort(eig(cov(scores)),'descend');
        totalvar(i)=sum(eig(cov(scores)));
        average_dist(i)=(mean((max(abs(fr_mat),[],1))));
        explained_mat(i,:)=explaineds(1:numPCs);
    end
    
    figure;
    co=colormap(jet(num_conds));
    set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
    plot(cumsum(explained_mat'),'linewidth',lw);
    legend(conds{1:num_conds});
    title('Variance Explained');
    set(gca,'fontsize',14);
    ylim([min(min(cumsum(explained_mat')))*0.95 100]);
    if ~isempty(totalvar)
        axes('Position',[.45 .25 .3 .3])
        box on
        plot(1:num_conds,totalvar,'-k','linewidth',lw);
        ax=gca;
        set(ax,'XTick',1:num_conds);
        set(ax,'XTickLabel',{'CC','Alt1','Alt2'});
        xlabel('Condition'); ylabel('Total Variance');set(gca,'fontsize',14);
        ylim([min(totalvar)*0.95 max(totalvar)*1.05]);
    end
    
    % Analyze variance explained by PCs of different conditions, with PC breakdown
    explaineds=cell(num_conds,1);
    condfig=nan(num_conds,1);
    for i=1:num_conds
        explaineds{i}=zeros(num_conds,numPCs);
        % Project each forward condition into alt and cc PCs
        condfig(i)=figure;
        co=colormap(lines(num_conds));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        PCs_cond=PCs{i};
        for c=1:num_conds
            if c==1
                fr_mat=cc_condense;
            elseif c==2
                fr_mat=a1_condense;
            elseif c==3
                fr_mat=a2_condense;
            end
            for cc=1:numPCs
                scores_thiscond=fr_mat*PCs_cond(:,1:cc);
                eigs_thiscond=eig(cov(scores_thiscond));
                explaineds{i}(c,cc)=sum(eigs_thiscond)/totalvar(c);
            end
            figure(condfig(i));
            if numPCs>2
                plot3(scores_thiscond(:,1),scores_thiscond(:,2),scores_thiscond(:,3),'linewidth',lw); hold on;
            elseif numPCs>1
                plot(scores_thiscond(:,1),scores_thiscond(:,2),'linewidth',lw); hold on;
            end
            fprintf('Percent variance explained by %s for %s : %d \n',conds{i},conds{c},round(100*explaineds{i}(c)))
        end
        figure(condfig(i));
        set(gca,'fontsize',18);
        title(strcat('All Cons proj in ',{' '},conds{i},{' '},' PC space'));
        xlabel(strcat(conds{i},{' '},' PC1')); ylabel(strcat(conds{i},{' '},' PC2')); zlabel(strcat(conds{i},{' '},' PC3'));
        legend(conds);
    end
    
    explaineds_6dim_condense=explaineds; % explaineds is cumulative
    
    %make plottable ver of cell array
    explaineds_2=explaineds;
    for i=1:length(explaineds)
        explaineds_2{i}=explaineds{i}(:,end);
    end
    
    PCsVarExplained=[{'CC';'Alt1';'Alt2'} explaineds_2(1:3)];
    PCsVarExplained_vals=PCsVarExplained(:,2);
    
    % PCA for all directions together, excluding INACT
    allA=[a1_condense; a2_condense; cc_condense; cc_condense];
    
    meanA=mean(allA,1);
    allA = allA-repmat(meanA,size(allA,1),1);
    [PCs,scores] = pca(allA);
    figure;
    co=colormap(lines(length(conds)));
    set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
    for i=1:length(conds)-1
        
        if i==1
            fr_mat=cc_condense;
        elseif i==2
            fr_mat=a1_condense;
        elseif i==3
            fr_mat=a2_condense;
        elseif i==4
            fr_mat=input(i).fr-repmat(mean(input(i).fr),size(input(i).fr,1),1);
        end
        
        scores=fr_mat*PCs(:,1:numPCs);
        if numPCs>2
            plot3(scores(:,1),scores(:,2),scores(:,3),'linewidth',lw); hold on;
        elseif numPCs>1
            plot(scores(:,1),scores(:,2),'linewidth',lw); hold on;
        end
    end
    set(gca,'fontsize',18);
    legend(conds);
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    grid on;
    %rotate3d on;
    
    % Make ratios (ratios for EMG are further below, after EMG PCA)
    neur_cc_ratio_across_6dims=nan*ones(1,6);
    neur_A1_ratio_across_6dims=nan*ones(1,6);
    
    for dim_num=1:6
        var_cc_in_A1 = explaineds_6dim_condense{2}(1,dim_num);
        var_cc_in_cc = explaineds_6dim_condense{1}(1,dim_num);
        neur_cc_ratio_across_6dims(dim_num)= var_cc_in_A1/var_cc_in_cc;
        
        var_A1_in_cc = explaineds_6dim_condense{1}(2,dim_num);
        var_A1_in_A1 = explaineds_6dim_condense{2}(2,dim_num);
        neur_A1_ratio_across_6dims(dim_num)= var_A1_in_cc/var_A1_in_A1;
    end
    
    
    %% run PCA on EMG
    
    warning('off','stats:pca:ColRankDefX')
    totalvar=zeros(length(conds),1);
    average_dist=zeros(length(conds),1);
    projs=cell(length(conds),1);
    PCs=cell(length(conds),1);
    explained_mat=zeros(length(conds),numPCs_emg);
    
    for i=1:length(conds)
        % this mean centers entire muscle matrix
        fr_mat=Data(i).emg-repmat(mean_EMG_across_conds,size(Data(i).emg,1),1);
        
        [PC,scores,latents,coeffs,explaineds_emg]=pca(fr_mat);
        PCs{i}=PC;
        projs{i}=scores(:,1:numPCs_emg);
        eigs=sort(eig(cov(scores)),'descend');
        totalvar(i)=sum(eig(cov(scores)));
        average_dist(i)=(mean((max(abs(fr_mat),[],1))));
        explained_mat(i,:)=explaineds_emg(1:numPCs_emg);
    end
    figure;
    co=colormap(jet(length(conds)));
    set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
    plot(cumsum(explained_mat'),'linewidth',lw);
    legend(conds{1:4});
    title('Variance Explained, EMG'); set(gca,'fontsize',14);
    ylim([min(min(cumsum(explained_mat')))*0.95 100]);
    if ~isempty(totalvar)
        axes('Position',[.45 .25 .3 .3])
        box on
        plot(1:length(conds),totalvar,'-k','linewidth',lw);
        ax=gca;
        set(ax,'XTick',1:length(conds));
        set(ax,'XTickLabel',{'CC','Alt1','Alt2','Inact'});
        xlabel('Condition'); ylabel('Total Variance');set(gca,'fontsize',14);
        ylim([min(totalvar)*0.95 max(totalvar)*1.05]);
    end
    
    
    %% Analyze variance explained by PCs of different conditions
    explaineds_emg=cell(length(conds),1);
    condfig=nan(length(conds),1);
    for i=1:length(conds)
        explaineds_emg{i}=zeros(length(conds),numPCs_emg);
        % Project each forward condition into alt and cc PCs
        condfig(i)=figure;
        co=colormap(lines(length(conds)));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        PCs_cond=PCs{i};
        for c=1:length(conds)
            fr_mat=Data(c).emg-repmat(mean_EMG_across_conds,size(Data(c).emg,1),1);
            
            for cc=1:numPCs_emg
                scores_thiscond=fr_mat*PCs_cond(:,1:cc);
                eigs_thiscond=eig(cov(scores_thiscond));
                explaineds_emg{i}(c,cc)=sum(eigs_thiscond)/totalvar(c);
            end
            
            figure(condfig(i));
            if numPCs_emg>2
                plot3(scores_thiscond(:,1),scores_thiscond(:,2),scores_thiscond(:,3),'linewidth',lw); hold on;
            elseif numPCs_emg>1
                plot(scores_thiscond(:,1),scores_thiscond(:,2),'linewidth',lw); hold on;
            end
            fprintf('Percent variance explained by %s for %s : %d \n',conds{i},conds{c},round(100*explaineds_emg{i}(c)),'EMG')
        end
        figure(condfig(i)); set(gca,'fontsize',18);
        title(strcat('All Conditions projected in ',{' '},conds{i},{' '},' PC space, EMG'));
        xlabel(strcat(conds{i},{' '},' PC1')); ylabel(strcat(conds{i},{' '},' PC2')); zlabel(strcat(conds{i},{' '},' PC3'));
        legend(conds{1:4});
    end
    
    %     figure;
    %     bar(1:length(conds),horzcat(explaineds_emg{:}));
    %     set(gca,'fontsize',18);legend(strcat(conds{1:4},' PCs'));
    %     ylim([0 1]);
    %     xlabel('Conditions');
    %     %xticklabels(conds);
    %     ax=gca; set(ax,'XTickLabel',conds);
    %     title('Variance explained in Different Conditions')
    
    PCsVarExplained=[{'CC';'Alt1';'Alt2';'Inact'} explaineds_emg];
    
    %% PCA for all directions together, INCLUDING INACT
    
    % get mean across all conditions
    allA=vertcat(Data.emg);
    meanA=mean(allA,1);
    
    %just use CC and A1 to define PC space
    allA=vertcat(Data(1:2).emg);
    % mean center using the mean across all conditions
    allA = allA-repmat(meanA,size(allA,1),1);
    
    [PCs,scores] = pca(allA);
    figure;
    co=colormap(lines(length(conds)));
    set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
    for i=1:length(conds)
        
        fr_mat=Data(i).emg-repmat(meanA,size(Data(i).emg,1),1);
        
        scores=fr_mat*PCs(:,1:numPCs_emg);
        if numPCs_emg>2
            plot3(scores(:,1),scores(:,2),scores(:,3),'linewidth',lw); hold on;
        elseif numPCs_emg>1
            plot(scores(:,1),scores(:,2),'linewidth',lw); hold on;
        end
    end
    set(gca,'fontsize',18); title('All Conditions projected in Alt1-CC PC space, EMG');
    legend(conds{1:4})
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    grid on;
    axis square;
    
    %% ratio variance explained for emg
    
    emg_cc_ratio_across_6dims=nan*ones(1,6);
    emg_A1_ratio_across_6dims=nan*ones(1,6);
    
    for dim_num=1:6
        var_cc_in_A1 = explaineds_emg{2}(1,dim_num);
        var_cc_in_cc = explaineds_emg{1}(1,dim_num);
        emg_cc_ratio_across_6dims(dim_num)= var_cc_in_A1/var_cc_in_cc;
        
        var_A1_in_cc = explaineds_emg{1}(2,dim_num);
        var_A1_in_A1 = explaineds_emg{2}(2,dim_num);
        emg_A1_ratio_across_6dims(dim_num)= var_A1_in_cc/var_A1_in_A1;
    end
    
    % plot the ratios
    figure;
    subplot(1,2,1); hold on;

    plot(1:6,emg_cc_ratio_across_6dims,'g');
    plot(1:6,neur_cc_ratio_across_6dims,'k');
    ylabel('Var of CC in top Alt1 PCs / Var of CC in top CC PCs');
    ylim([0 1]);
    title('CC');
    
    subplot(1,2,2); hold on;
    plot(1:6,emg_A1_ratio_across_6dims,'g');
    plot(1:6,neur_A1_ratio_across_6dims,'k');
    legend({'EMG','neural'});
    ylabel('Var of Alt1 in top CC PCs / Var of Alt1 in top Alt1 PCs');
    ylim([0 1]);
    title('Alt1');

    %% similarity index on EMG

    align_tags={'A1 and CC','A2 and CC','A1 and A2','A1 and Inact','CC and Inact'};
    
    input=Data;
    PCAdim_emg=6;
    align = zeros(5,1);
    cc_mat=input(1).emg-repmat(mean(input(1).emg),size(input(1).emg,1),1);
    a1_mat=input(2).emg-repmat(mean(input(2).emg),size(input(2).emg,1),1);
    a2_mat=input(3).emg-repmat(mean(input(3).emg),size(input(3).emg,1),1);
    inact_mat=input(4).emg-repmat(mean(input(4).emg),size(input(4).emg,1),1);
    
    [PCsCC,~,evaluesCC] = pca(cc_mat);
    [PCsA1,~,evaluesA1] = pca(a1_mat);
    [PCsA2,~,evaluesA2] = pca(a2_mat);
    [PCsInact,~,evaluesInact] = pca(inact_mat);
    
    CCinA1 = (trace(PCsA1(:,1:PCAdim_emg)'*cov(cc_mat)*PCsA1(:,1:PCAdim_emg))) / (trace(PCsCC(:,1:PCAdim_emg)'*cov(cc_mat)*PCsCC(:,1:PCAdim_emg)));
    A1inCC = (trace(PCsCC(:,1:PCAdim_emg)'*cov(a1_mat)*PCsCC(:,1:PCAdim_emg))) / (trace(PCsA1(:,1:PCAdim_emg)'*cov(a1_mat)*PCsA1(:,1:PCAdim_emg)));
    align(1) = mean([CCinA1 A1inCC]);
    
    CCinA2 = (trace(PCsA2(:,1:PCAdim_emg)'*cov(cc_mat)*PCsA2(:,1:PCAdim_emg))) / (trace(PCsCC(:,1:PCAdim_emg)'*cov(cc_mat)*PCsCC(:,1:PCAdim_emg)));
    A2inCC = (trace(PCsCC(:,1:PCAdim_emg)'*cov(a2_mat)*PCsCC(:,1:PCAdim_emg))) / (trace(PCsA2(:,1:PCAdim_emg)'*cov(a1_mat)*PCsA2(:,1:PCAdim_emg)));
    align(2) = mean([CCinA2 A2inCC]);
    
    A2inA1 = (trace(PCsA1(:,1:PCAdim_emg)'*cov(a2_mat)*PCsA1(:,1:PCAdim_emg))) / (trace(PCsA2(:,1:PCAdim_emg)'*cov(a2_mat)*PCsA2(:,1:PCAdim_emg)));
    A1inA2 = (trace(PCsA2(:,1:PCAdim_emg)'*cov(a1_mat)*PCsA2(:,1:PCAdim_emg))) / (trace(PCsA1(:,1:PCAdim_emg)'*cov(a1_mat)*PCsA1(:,1:PCAdim_emg)));
    align(3) = mean([A2inA1 A1inA2]);
    
    InactinA1 = (trace(PCsA1(:,1:PCAdim_emg)'*cov(inact_mat)*PCsA1(:,1:PCAdim_emg))) / (trace(PCsInact(:,1:PCAdim_emg)'*cov(inact_mat)*PCsInact(:,1:PCAdim_emg)));
    A1inInact = (trace(PCsInact(:,1:PCAdim_emg)'*cov(a1_mat)*PCsInact(:,1:PCAdim_emg))) / (trace(PCsA1(:,1:PCAdim_emg)'*cov(a1_mat)*PCsA1(:,1:PCAdim_emg)));
    align(4) = mean([InactinA1 A1inInact]);
    
    InactinCC = (trace(PCsCC(:,1:PCAdim_emg)'*cov(inact_mat)*PCsCC(:,1:PCAdim_emg))) / (trace(PCsInact(:,1:PCAdim_emg)'*cov(inact_mat)*PCsInact(:,1:PCAdim_emg)));
    CCinInact = (trace(PCsInact(:,1:PCAdim_emg)'*cov(cc_mat)*PCsInact(:,1:PCAdim_emg))) / (trace(PCsCC(:,1:PCAdim_emg)'*cov(cc_mat)*PCsCC(:,1:PCAdim_emg)));
    align(5) = mean([InactinA1 CCinInact]);
    
    figure; hold on;
    plot([1:5],align,'*k');
    
    ylim([0 1]);
    xlim([0.5 5.5]);
    xticks([1:5])
    xticklabels(align_tags);
    ylabel('Alignment index');
    
    align_emg=align;
    title('Subspace alignment');
    
    %% similarity index on neural data

    align_tags={'A1 and CC','A2 and CC','A1 and A2','A1 and Inact','CC and Inact'};
    
    for ii=1:2
        if ii==1
            input=Data;
        elseif ii==2
            input=Data_pyr;
        end
        
        align = zeros(5,1);
        
        cc_mat=input(1).fr-repmat(mean(input(1).fr),size(input(1).fr,1),1);
        a1_mat=input(2).fr-repmat(mean(input(2).fr),size(input(2).fr,1),1);
        a2_mat=input(3).fr-repmat(mean(input(3).fr),size(input(3).fr,1),1);
        inact_mat=input(4).fr-repmat(mean(input(4).fr),size(input(4).fr,1),1);
        
        [PCsCC,~,evaluesCC] = pca(cc_mat);
        [PCsA1,~,evaluesA1] = pca(a1_mat);
        [PCsA2,~,evaluesA2] = pca(a2_mat);
        [PCsInact,~,evaluesInact] = pca(inact_mat);
        
        CCinA1 = (trace(PCsA1(:,1:PCAdim)'*cov(cc_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        A1inCC = (trace(PCsCC(:,1:PCAdim)'*cov(a1_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(1) = mean([CCinA1 A1inCC]);
        
        CCinA2 = (trace(PCsA2(:,1:PCAdim)'*cov(cc_mat)*PCsA2(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        A2inCC = (trace(PCsCC(:,1:PCAdim)'*cov(a2_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsA2(:,1:PCAdim)'*cov(a1_mat)*PCsA2(:,1:PCAdim)));
        align(2) = mean([CCinA2 A2inCC]);
        
        A2inA1 = (trace(PCsA1(:,1:PCAdim)'*cov(a2_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsA2(:,1:PCAdim)'*cov(a2_mat)*PCsA2(:,1:PCAdim)));
        A1inA2 = (trace(PCsA2(:,1:PCAdim)'*cov(a1_mat)*PCsA2(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(3) = mean([A2inA1 A1inA2]);
        
        InactinA1 = (trace(PCsA1(:,1:PCAdim)'*cov(inact_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsInact(:,1:PCAdim)'*cov(inact_mat)*PCsInact(:,1:PCAdim)));
        A1inInact = (trace(PCsInact(:,1:PCAdim)'*cov(a1_mat)*PCsInact(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(4) = mean([InactinA1 A1inInact]);
        
        InactinCC = (trace(PCsCC(:,1:PCAdim)'*cov(inact_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsInact(:,1:PCAdim)'*cov(inact_mat)*PCsInact(:,1:PCAdim)));
        CCinInact = (trace(PCsInact(:,1:PCAdim)'*cov(cc_mat)*PCsInact(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        align(5) = mean([InactinA1 CCinInact]);
        
        if ii==1
            align_all=align;
        elseif ii==2
            align_pyr=align;
        end
        
    end
    
    figure; hold on;
    plot([1:5],align_all,'*k');
    plot([1:5],align_pyr,'*b');
    ylim([0 1]);
    xlim([0.5 5.5]);
    xticks([1:5])
    xticklabels(align_tags);
    ylabel('Alignment index');
    title('Subspace alignment');
    
    %% RE-RUN ANALYSES ON MUCLE-RELATED NEURAL DATA
    % but NOT the Inact, keep that data "uncleaned" bc there's no muscle
    % activity to filter by
    
    % similarity index
    for ii=1:2
        if ii==1
            input=Data;
        elseif ii==2
            input=Data_pyr;
        end
        align = zeros(5,1);
        cc_mat=input(1).clean_data-repmat(mean(input(1).clean_data),size(input(1).clean_data,1),1);
        a1_mat=input(2).clean_data-repmat(mean(input(2).clean_data),size(input(2).clean_data,1),1);
        a2_mat=input(3).clean_data-repmat(mean(input(3).clean_data),size(input(3).clean_data,1),1);
        inact_mat=input(4).fr-repmat(mean(input(4).fr),size(input(4).fr,1),1); % UNCLEANED INACT
        
        [PCsCC,~,evaluesCC] = pca(cc_mat);
        [PCsA1,~,evaluesA1] = pca(a1_mat);
        [PCsA2,~,evaluesA2] = pca(a2_mat);
        [PCsInact,~,evaluesInact] = pca(inact_mat);
        
        CCinA1 = (trace(PCsA1(:,1:PCAdim)'*cov(cc_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        A1inCC = (trace(PCsCC(:,1:PCAdim)'*cov(a1_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(1) = mean([CCinA1 A1inCC]);
        
        CCinA2 = (trace(PCsA2(:,1:PCAdim)'*cov(cc_mat)*PCsA2(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        A2inCC = (trace(PCsCC(:,1:PCAdim)'*cov(a2_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsA2(:,1:PCAdim)'*cov(a1_mat)*PCsA2(:,1:PCAdim)));
        align(2) = mean([CCinA2 A2inCC]);
        
        A2inA1 = (trace(PCsA1(:,1:PCAdim)'*cov(a2_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsA2(:,1:PCAdim)'*cov(a2_mat)*PCsA2(:,1:PCAdim)));
        A1inA2 = (trace(PCsA2(:,1:PCAdim)'*cov(a1_mat)*PCsA2(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(3) = mean([A2inA1 A1inA2]);
        
        InactinA1 = (trace(PCsA1(:,1:PCAdim)'*cov(inact_mat)*PCsA1(:,1:PCAdim))) / (trace(PCsInact(:,1:PCAdim)'*cov(inact_mat)*PCsInact(:,1:PCAdim)));
        A1inInact = (trace(PCsInact(:,1:PCAdim)'*cov(a1_mat)*PCsInact(:,1:PCAdim))) / (trace(PCsA1(:,1:PCAdim)'*cov(a1_mat)*PCsA1(:,1:PCAdim)));
        align(4) = mean([InactinA1 A1inInact]);
        
        InactinCC = (trace(PCsCC(:,1:PCAdim)'*cov(inact_mat)*PCsCC(:,1:PCAdim))) / (trace(PCsInact(:,1:PCAdim)'*cov(inact_mat)*PCsInact(:,1:PCAdim)));
        CCinInact = (trace(PCsInact(:,1:PCAdim)'*cov(cc_mat)*PCsInact(:,1:PCAdim))) / (trace(PCsCC(:,1:PCAdim)'*cov(cc_mat)*PCsCC(:,1:PCAdim)));
        align(5) = mean([InactinA1 CCinInact]);
        
        if ii==1
            align_all_clean=align;
        elseif ii==2
            align_pyr_clean=align;
        end
        
    end
    
    figure; hold on;
    plot([1:5],align_all_clean,'*k');
    plot([1:5],align_pyr_clean,'*b');
    legend('all','pyr');
    ylim([0 1]);
    xlim([0.5 5.5]);
    xticks([1:5]);
    xticklabels(align_tags);
    ylabel('Alignment index');
    title(strcat(mouseName, ' subspace alignment cleaned'));
    
    
    %% PCA on muscle-related neural data
    for kk=1:2
        if kk==1
            input=Data;
        elseif kk==2
            input=Data_pyr;
        end
        warning('off','stats:pca:ColRankDefX')
        totalvar=zeros(length(conds),1);
        average_dist=zeros(length(conds),1);
        projs=cell(length(conds),1);
        PCs=cell(length(conds),1);
        explained_mat=zeros(length(conds),numPCs);
        
        %since using un-filtered Inact, only run this loop on 3 conds, 4th below
        for i=1:length(conds)-1
            fr_mat=input(i).clean_data-repmat(mean(input(i).clean_data),size(input(i).clean_data,1),1);
            [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
            PCs{i}=PC;
            projs{i}=scores(:,1:numPCs);
            eigs=sort(eig(cov(scores)),'descend');
            totalvar(i)=sum(eig(cov(scores)));
            average_dist(i)=(mean((max(abs(fr_mat),[],1))));
            explained_mat(i,:)=explaineds(1:numPCs);
        end
        
        % for un-filtered inact
        for i=4
            %fr_mat=input(i).clean_data-repmat(mean(input(i).clean_data),size(input(i).clean_data,1),1);
            fr_mat=input(i).fr-repmat(mean(input(i).fr),size(input(i).fr,1),1);
            [PC,scores,latents,coeffs,explaineds]=pca(fr_mat);
            PCs{i}=PC;
            projs{i}=scores(:,1:numPCs);
            eigs=sort(eig(cov(scores)),'descend');
            totalvar(i)=sum(eig(cov(scores)));
            average_dist(i)=(mean((max(abs(fr_mat),[],1))));
            explained_mat(i,:)=explaineds(1:numPCs);
        end
        
        figure;
        co=colormap(jet(length(conds)));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        plot(cumsum(explained_mat'),'linewidth',lw);
        legend(conds{1:4});
        if kk==1
            title('Variance Explained (Neural) ALL cleaned');
        elseif kk==2
            title('Variance Explained (Neural) PYR cleaned');
        end
        set(gca,'fontsize',14);
        ylim([min(min(cumsum(explained_mat')))*0.95 100]);
        if ~isempty(totalvar)
            axes('Position',[.45 .25 .3 .3])
            box on
            plot(1:length(conds),totalvar,'-k','linewidth',lw);
            ax=gca;
            set(ax,'XTick',1:length(conds));
            set(ax,'XTickLabel',{'CC','Alt1','Alt2','Inact'});
            xlabel('Condition'); ylabel('Total Variance');set(gca,'fontsize',14);
            ylim([min(totalvar)*0.95 max(totalvar)*1.05]);
        end
        
        % Analyze variance explained by PCs of different conditions, INCLUDING INACT WITH PC breakdown
        explaineds=cell(length(conds),1);
        condfig=nan(length(conds),1);
        for i=1:length(conds)
            explaineds{i}=zeros(length(conds),numPCs);
            % Project each forward condition into alt and cc PCs
            condfig(i)=figure;
            co=colormap(lines(length(conds)));
            set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
            PCs_cond=PCs{i};
            for c=1:length(conds)
                fr_mat=input(c).clean_data-repmat(mean(input(c).clean_data),size(input(c).clean_data,1),1);
                for cc=1:numPCs
                    scores_thiscond=fr_mat*PCs_cond(:,1:cc);
                    eigs_thiscond=eig(cov(scores_thiscond));
                    explaineds{i}(c,cc)=sum(eigs_thiscond)/totalvar(c);
                end
                figure(condfig(i));
                if numPCs>2
                    plot3(scores_thiscond(:,1),scores_thiscond(:,2),scores_thiscond(:,3),'linewidth',lw); hold on;
                elseif numPCs>1
                    plot(scores_thiscond(:,1),scores_thiscond(:,2),'linewidth',lw); hold on;
                end
                fprintf('Percent variance explained by %s for %s : %d \n',conds{i},conds{c},round(100*explaineds{i}(c)))
            end
            figure(condfig(i));
            set(gca,'fontsize',18);
            
            if kk==1
                title(strcat('All Cons proj in ',{' '},conds{i},{' '},' PC space ALL cleaned'));
            elseif kk==2
                title(strcat('All Cons proj in ',{' '},conds{i},{' '},' PC space PYR cleaned'));
            end
            
            xlabel(strcat(conds{i},{' '},' PC1')); ylabel(strcat(conds{i},{' '},' PC2')); zlabel(strcat(conds{i},{' '},' PC3'));
            legend(conds);
            grid on;
            axis square;
        end
        if kk==1
            explaineds_all_clean=explaineds;
        elseif kk==2
            explaineds_pyr_clean=explaineds;
        end
        
        %            %make plottable ver of cell array
        %            explaineds_2=explaineds;
        %            for i=1:length(explaineds)
        %                explaineds_2{i}=explaineds{i}(:,end);
        %            end
        %
        %            figure;
        %            bar(1:length(conds),horzcat(explaineds_2{:})');
        %            set(gca,'fontsize',18);
        %            legend(strcat(conds,' PCs'));
        %            ylim([0 1]);
        %            xlabel('Conditions');
        %            %xticklabels(conds);
        %            ax=gca; set(ax,'XTickLabel',conds);
        %            if kk==1
        %                title('Var explained in Diff Conds ALL cleaned');
        %            elseif kk==2
        %                title('Var explained in Diff Conds PYR cleaned');
        %            end
        
        figure;
        subplot(1,4,1);
        plot(explaineds{1}');
        ylim([0 1]);
        title('CC PCs');
        subplot(1,4,2);
        plot(explaineds{2}');
        ylim([0 1]);
        title('A1 PCs');
        subplot(1,4,3);
        plot(explaineds{3}');
        ylim([0 1]);
        title('A2 PCs');
        subplot(1,4,4);
        plot(explaineds{4}');
        ylim([0 1]);
        legend('cc','a1','a2','inact');
        title('Inact PCs');
        if kk==1
            suptitle('ALL cleaned');
        elseif kk==2
            suptitle('PYRs cleaned');
        end
        
        % PCA for all directions together, INCLUDING INACT
        allA=vertcat(input([1:2]).clean_data); %make PC space w/o Alt2 4 conds
        meanA=mean(allA,1);
        allA = allA-repmat(meanA,size(allA,1),1);
        [PCs,scores] = pca(allA);
        figure;
        co=colormap(lines(length(conds)));
        set(gca, 'ColorOrder', co, 'NextPlot', 'replacechildren');
        for i=1:length(conds)
            if size(meanA)==size(input(i).clean_data)
                fr_mat=input(i).clean_data-meanA;
            else
                fr_mat=input(i).clean_data-repmat(meanA,size(input(i).clean_data,1),1);
            end
            scores=fr_mat*PCs(:,1:numPCs);
            if numPCs>2
                plot3(scores(:,1),scores(:,2),scores(:,3),'linewidth',lw); hold on;
            elseif numPCs>1
                plot(scores(:,1),scores(:,2),'linewidth',lw); hold on;
            end
        end
        set(gca,'fontsize',18);
        if kk==1
            title('All Conds proj in PC space ALL cleaned');
            if mouse==1
                xlim([-60 60]);
                ylim([-60 60]);
                zlim([-60 60]);
                xticks([-60 0 60])
                yticks([-60 0 60])
                zticks([-60 0 60])
            end
        elseif kk==2
            title('All Conds proj in PC space PYR cleaned');
            if mouse==1
                xlim([-60 60]);
                ylim([-60 60]);
                zlim([-60 60]);
                xticks([-60 0 60])
                yticks([-60 0 60])
                zticks([-60 0 60])
            end
        end
        legend(conds);
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        grid on;
        axis square;
    end
    
    %% save analysis data
    %         Data_fr_and_clean_all=struct;
    %         Data_fr_and_clean_pyr=struct;
    
    for i=1:4
        Data_fr_and_clean_all(i).fr=Data(i).fr;
        Data_fr_and_clean_all(i).emg=Data(i).emg;
        Data_fr_and_clean_all(i).clean_data=Data(i).clean_data;
        Data_fr_and_clean_pyr(i).fr=Data_pyr(i).fr;
        Data_fr_and_clean_pyr(i).emg=Data_pyr(i).emg;
        Data_fr_and_clean_pyr(i).clean_data=Data_pyr(i).clean_data;
    end
    
    filename=strcat(mouse_name,'_analysis_data_final.mat');
    
    if saveon
        save(filename,'mouse_data','conds','muscle_names','neurons_total','neurons_select','fracfrFactor','minfr','slopeThresh',...
            'num_perms_lambda','numPCs','numPCs_emg','numPCs_forClean','trial_idxs','neuron_origin_ID','drop_sessions','totalfrCell',...
            'slopes_data','width_IN','width_pyr','Data_fr_and_clean_all','Data_fr_and_clean_pyr',...
            'Data_subtypes','width_15ms_keep','waveforms_select',...
            'explaineds_all','explaineds_pyr','explaineds_emg','explaineds_all_clean','explaineds_pyr_clean',...
            'align_tags','align_all','align_pyr','align_all_clean','align_pyr_clean','align_emg',...
            'ortho_control_neuro','ortho_control_emg',...
            'neur_cc_ratio_across_6dims','neur_A1_ratio_across_6dims',...
            'emg_cc_ratio_across_6dims','emg_A1_ratio_across_6dims',...
            'R2_conds_joint_scale_all','R2_conds_joint_scale_pyr',...
            'Y_real_CCmodel_all','Y_fit_CCmodel_all','M_CCmodel_all','Test_R2_CCmodel_all','lambda_opt_CCmodel_all',...
            'Y_real_A1model_all','Y_fit_A1model_all','M_A1model_all','Test_R2_A1model_all','lambda_opt_A1model_all',...
            'Y_real_CCmodel_pyr','Y_fit_CCmodel_pyr','M_CCmodel_pyr','Test_R2_CCmodel_pyr','lambda_opt_CCmodel_pyr',...
            'Y_real_A1model_pyr','Y_fit_A1model_pyr','M_A1model_pyr','Test_R2_A1model_pyr','lambda_opt_A1model_pyr',...
            'M_A2model_all','M_A2model_pyr',...
            'Y_fit_CCmodel_Xperms_all','Y_fit_A1model_Xperms_all','Y_fit_CCmodel_Xperms_pyr','Y_fit_A1model_Xperms_pyr',...
            'Y_real_joint_all','Y_fit_joint_all','M_joint_all','Test_R2_joint_all','lambda_opt_joint_all',...
            'Y_real_joint_pyr','Y_fit_joint_pyr','M_joint_pyr','Test_R2_joint_pyr','lambda_opt_joint_pyr',...
            'dropNeurs','Test_R2_ortho_cc','Test_R2_ortho_cc_pyr','-v7.3');
    end
    
    clearvars -except mouse
end
