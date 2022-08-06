%%% Overlay the good ROIs onto static ims to classify cells as rabies-positive or not
close all
clear all

%parameters
orig_line_leng = 256; %In pixels
frames2avg = 128;
addl_corr = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%Initialize
next_fig_num = 1;

%Define files to use
filename(1).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D19';
filename(2).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D19';
filename(3).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D20';
filename(4).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D49';
filename(5).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D49';
filename(6).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D49';
filename(7).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D49';
filename(8).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D50';
filename(9).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D50';
filename(10).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D50';
filename(11).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D50';
filename(12).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D50';
% filename(13).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D28';
% filename(14).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D28';
% filename(15).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D28';
% filename(16).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D28';
% filename(17).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D28';
% filename(18).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D26';

filename(1).mov = 'TSeries-11012017-D19-001_mov_nonrigid_corr_tot.mat';
filename(2).mov = 'TSeries-11022017-D19-001_mov_nonrigid_corr_tot.mat';
filename(3).mov = 'TSeries-11032017-D20-001_mov_nonrigid_corr_tot.mat';
filename(4).mov = 'TSeries-04192018-D49-001_mov_nonrigid_corr_tot.mat';
filename(5).mov = 'TSeries-04222018-D49-001_mov_nonrigid_corr_tot.mat';
filename(6).mov = 'TSeries-04232018-D49-001_mov_nonrigid_corr_tot.mat';
filename(7).mov = 'TSeries-04242018-D49-001_mov_nonrigid_corr_tot.mat';
filename(8).mov = 'TSeries-04192018-D50-001_mov_nonrigid_corr_tot.mat';
filename(9).mov = 'TSeries-04212018-D50-001_mov_nonrigid_corr_tot.mat';
filename(10).mov = 'TSeries-04222018-D50-001_mov_nonrigid_corr_tot.mat';
filename(11).mov = 'TSeries-04232018-D50-001_mov_nonrigid_corr_tot.mat';
filename(12).mov = 'TSeries-04252018-D50-001_mov_nonrigid_corr_tot.mat';
% filename(13).mov = 'TSeries-11282017-D28-001_mov_nonrigid_corr_tot.mat';
% filename(14).mov = 'TSeries-11302017-D28-001_mov_nonrigid_corr_tot.mat';
% filename(15).mov = 'TSeries-12012017-D28-001_mov_nonrigid_corr_tot.mat';
% filename(16).mov = 'TSeries-12022017-D28-001_mov_nonrigid_corr_tot.mat';
% filename(17).mov = 'TSeries-12022017-D28-002_mov_nonrigid_corr_tot.mat';
% filename(18).mov = 'TSeries-11282017-D26-001_mov_nonrigid_corr_tot.mat';

filename(1).cnmf = 'TSeries-11012017-D19-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(2).cnmf = 'TSeries-11022017-D19-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(3).cnmf = 'TSeries-11032017-D20-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(4).cnmf = 'TSeries-04192018-D49-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(5).cnmf = 'TSeries-04222018-D49-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(6).cnmf = 'TSeries-04232018-D49-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(7).cnmf = 'TSeries-04242018-D49-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(8).cnmf = 'TSeries-04192018-D50-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(9).cnmf = 'TSeries-04212018-D50-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(10).cnmf = 'TSeries-04222018-D50-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(11).cnmf = 'TSeries-04232018-D50-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(12).cnmf = 'TSeries-04252018-D50-001_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(13).cnmf = 'TSeries-11282017-D28-001_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(14).cnmf = 'TSeries-11302017-D28-001_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(15).cnmf = 'TSeries-12012017-D28-001_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(16).cnmf = 'TSeries-12022017-D28-001_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(17).cnmf = 'TSeries-12022017-D28-002_CNMF_results_1nb_1_5_p1_fin.mat';
% filename(18).cnmf = 'TSeries-11282017-D26-001_CNMF_results_1nb_1_5_p1_fin.mat';

%Main FOR loop
for i = 1:12
    i
    names = strsplit(filename(i).mov,'-');
    names2 = strsplit(names{4},'_');
    savename = [char(names{3}) '_' char(names{2}) '_' char(names2{1}) '_mean_ims_' num2str(frames2avg)];
    savename_Ybi = [char(names{3}) '_' char(names{2}) '_' char(names2{1}) '_Y_bi'];
    
    %Load data
    load_time = tic;
    cd(filename(i).ts_dir)
    load(filename(i).mov)
    
%     if i == 9
%         Y(:,:,1:18000) = mov_nonrigid_corr(:,:,1:18000);
%         Y(:,:,18001:36000) = mov_nonrigid_corr(:,:,36001:54000);
%     else
        Y = mov_nonrigid_corr;
%     end
    clear mov_nonrigid_corr
    
    %Convert to single
    if ~isa(Y,'single');   
        Y = single(Y);
    end
    load_time = toc(load_time)
    
    [~,Y_bi] = correct_bidirectional_offset(Y,size(Y,3));
    clear Y
    
    if addl_corr(i) == 1
        options_rigid = NoRMCorreSetParms('d1',size(Y_bi,1),'d2',size(Y_bi,2),'bin_width',50,'max_shift',15,'us_fac',50);
        mot_corr_time = tic;
        [Y_bi,~,~] = normcorre(Y_bi,options_rigid);
        toc(mot_corr_time)
    end
    
    % Trim off edges with incomplete sampling
    load(filename(i).cnmf,'rows','cols');
    Y_bi(rows,:,:) = [];
    Y_bi(:,cols,:) = [];
    
    %Make Mean Images
    num_frames = size(Y_bi,3);
    mean_im_start = mean(Y_bi(:,:,1:frames2avg),3);
    mean_im_mid = mean(Y_bi(:,:,ceil(num_frames/2)+1:ceil(num_frames/2)+frames2avg),3);
    mean_im_end = mean(Y_bi(:,:,end-frames2avg:end),3);
    
    save(savename,'mean_im_start','mean_im_mid','mean_im_end');
%     save(savename_Ybi,'Y_bi','-v7.3')
    clear Y_bi
end

            
            

