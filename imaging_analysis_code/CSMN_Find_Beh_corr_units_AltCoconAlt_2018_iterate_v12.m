%%%% Process CSV files with EMG, etc. and output from Simons processing
%%%According to Bruker, collection begins with first galvo rotation and
%%%does not skip any, such that only galvo rotations on the end of the time
%%%series should be ignored. Verify with Pockels blanking signal in the
%%%future
close all
clear
runtime = tic;

animal = 'D16';
data_path = 'C:\Users\jam5064\Desktop\Current Data\D16_imaging_data\TSeries\TSeries-11222017-PM';
csv_filenames = {'TSeries-11222017-2205-D16-001_Cycle00001_VoltageRecording_001.csv','TSeries-11222017-2205-D16-003_Cycle00001_VoltageRecording_001.csv'};
Ca_filename = 'TSeries-11222017b-D16-001_CNMF_results_1nb_1_5_p1.mat';
mph = 0.3; %MinPeakHeight

names = strsplit(csv_filenames{1},'-');
voltages_filename = [animal '_' char(names{2}) '_voltages'];
segments_filename = [animal '_' char(names{2}) '_segments'];
names = strsplit(Ca_filename,'-');
names2 = strsplit(char(names{4}),'.');
align_filename = [animal '_' char(names{2}) '_align_v12_' char(names2{1})];

%Initialize
next_fig_num = 1;

cd(data_path)
load(voltages_filename);
load(segments_filename)
load(Ca_filename);
num_cells = size(A,2);

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
num_stds = 3;
noise_comp_cuttoff = round(num_cells/2);
align_pc_comp = zeros(2,3,4,3);
good_cells_full = [];
var_exp = cell(2,4,3);
for PCAdim = [5 6 7]
    %calcium - good cells define with spiking tho
    coeff_sp_alt1 = pca(alt1_means(:,:,3));
    noise_comps = abs(coeff_sp_alt1(:,PCAdim+1:end));
    cell_comps_mean = mean(noise_comps,2);
    cell_comps_std = std(noise_comps,0,2);
    comp_alt = repmat(cell_comps_mean + num_stds*cell_comps_std,1,PCAdim);
    
    coeff_sp_cocon = pca(cocon_means(:,:,3));
    noise_comps = abs(coeff_sp_cocon(:,PCAdim+1:end));
    cell_comps_mean = mean(noise_comps,2);
    cell_comps_std = std(noise_comps,0,2);
    comp_cocon = repmat(cell_comps_mean + num_stds*cell_comps_std,1,PCAdim);
    
    good_cells = find(sum(abs(coeff_sp_alt1(:,1:PCAdim))>comp_alt,2) + sum(abs(coeff_sp_cocon(:,1:PCAdim))>comp_cocon,2));
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
