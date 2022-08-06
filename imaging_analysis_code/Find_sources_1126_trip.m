clear; 
close all;
total_time = tic;

%% Set params and load data
save_results = 1;
play_movies = 1;
addl_corr = 0;
data_path = ('C:\Users\jam5064\Desktop\Current Data\D50');
filename = 'TSeries-04252018-D50-001_mov_nonrigid_corr_tot.mat'
names = strsplit(filename,'_');
savename_CNMF_results = [char(names{1}) '_CNMF_results_1nb_1_5_p1_fin']; %using only 1 bg comp may make a big diff??
frames2avg = 100;
next_fig_num = 1;

%%% Set CNMF parameters
K = 200;    % number of components to be found
tau = 1;     % std of gaussian kernel (size of neuron) - must be integer for manual refinement
p = 1; 

gcp;  % start cluster
addpath(genpath('utilities'));
addpath(genpath('deconvolution'));

%%% Load Data
load_time = tic;

cd(data_path)
load(filename)
Y = mov_nonrigid_corr;
clear mov_nonrigid_corr

if ~isa(Y,'single');   % convert to single
    Y = single(Y);
end

load_time = toc(load_time)


%% Main

%%% Start run timer
corr_time = tic;

%correct for bidirectional
[~,Y_bi] = correct_bidirectional_offset(Y,size(Y,3));
figure(next_fig_num); next_fig_num = next_fig_num+1;
imagesc(mean(Y_bi,3)); axis image %sanity check
title('Corrected bidirectional')
pause(5)

% make meanMovIm for alignment before trimming (I think this is nec)
image_data_struct.meanMovIm = mean(Y_bi(:,:,end-(frames2avg-1):end),3);

%play back as sanity check
if play_movies
    pause
    q = quantile(Y_bi(:),[0.005 0.995]);
    figure(next_fig_num); next_fig_num = next_fig_num+1;
    play_movie_array(Y_bi,q(1),q(2),0.001)
    close all
	addl_corr = input('Residual motion?')
end

if addl_corr
    options_rigid = NoRMCorreSetParms('d1',size(Y_bi,1),'d2',size(Y_bi,2),'bin_width',50,'max_shift',15,'us_fac',50);
    mot_corr_time = tic; 
    [Y_bi,~,~] = normcorre(Y_bi,options_rigid); 
    toc(mot_corr_time)
    if play_movies
        pause
        q = quantile(Y_bi(:),[0.005 0.995]);
        play_movie_array(Y_bi,q(1),q(2))
    end
end

% Trim off edges with incomplete sampling
rows = [];
cols = [];
for i = 10:10:size(Y_bi,3)
    rows = unique([rows; find(sum(sum(Y_bi(:,:,i-9:i),3),2)==0)]);
    cols = unique([cols find(sum(sum(Y_bi(:,:,i-9:i),3))==0)]);
end
rows
cols
Y_bi(rows,:,:) = [];
Y_bi(:,cols,:) = [];

figure(next_fig_num); next_fig_num = next_fig_num+1;
imagesc(mean(Y_bi,3)); 
title('Corrected after trimming')
axis image 
pause(5)

correction_time = toc(corr_time)/60

init_time = tic; 

% Calc size
[d1,d2,T] = size(Y_bi);
d = d1*d2; % total number of pixels

%%% update spatial components
options = CNMFSetParms(...
    'd1',d1,'d2',d2,...                         % dimensionality of the FOV
    'tsub',4,...
    'bSiz',1,...                                % expand kernel for HALS growing (default: 3)
    'nb',1,.....                                % number of background components (default: 1) %using only 1 bg comp may make a big diff??
    'maxthr',0.5,...                            % threshold of max value below which values are discarded (default: 0.25)
    'p',p,...                                   % order of AR dynamics
    'gSig',1,...                                % half size of neuron
    'medw',[2 2],...                            % size of median filter (default: [3,3])
    'spatial_method','constrained',...
    'thr_method','nrg',...
    'nrgthr',0.95,...
    'search_method','ellipse',...
    'min_size',1,...
    'max_size',5,...
    'deconv_method','MCMC' ...
    );

%%% Data pre-processing
%%% Find initial ROIs with iterative thresholding
[P,Y_bi] = preprocess_data(Y_bi,p);
use_mm = image_data_struct.meanMovIm;
use_mm(rows,:,:) = [];
use_mm(:,cols,:) = [];
[miri_masks, P.ROI_list] = get_dendrite_rois_onepass(use_mm);
P.ROI_list = [P.ROI_list(:,2) P.ROI_list(:,1)]; %need to reverse

[Ain,Cin,bin,fin,center] = initialize_components(Y_bi,K,tau,options,P);

%%% Want to refine components to break apart large chunks
[Ain,Cin,center] = manually_refine_components(Y_bi,Ain,Cin,center,use_mm,tau,options);

%%% display centers of found components
figure(next_fig_num); next_fig_num = next_fig_num+1;
plot_contours(Ain,use_mm,options,0); % contour plot of spatial footprints
title(['After initialize, i = ' num2str(i) ', j = ' num2str(j)])

initialize_time = toc(init_time)/60

spat_time = tic; 

Yr = reshape(Y_bi,d,T);
[A,b,C] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
R = (double(Yr)'*A)'; % raw roi time series

figure(next_fig_num); next_fig_num = next_fig_num+1;
plot_contours(A,use_mm,options,0); % contour plot of spatial footprints
title(['After spatial update, i = ' num2str(i) ', j = ' num2str(j)])

spatial_time = toc(spat_time)/60

bl_time = tic;

%fix baseline
Rsub = R;
num_cells = size(R,1);
num_samples = size(R,2);
bin_size = 25;
wb=waitbar(0,'baseline correction');
for i = 1:num_cells
    for j = 1:bin_size:num_samples-bin_size+1
        Rsub(i,j:j+bin_size-1) = Rsub(i,j:j+bin_size-1) - quantile(Rsub(i,j:j+bin_size-1),0.1);
        waitbar(i/num_cells,wb,'baseline correction');
    end
end
close(wb)
baseline_time = toc(bl_time)/60

temp_time = tic;

%deconv
if isempty(gcp)
    parpool;
end
Rdcv = Rsub*0;
S = Rsub*0;
params.B = 300;
params.Nsamples = 400;
params.p = P.p;
params.bas_nonneg = options.bas_nonneg;
parfor i = 1:num_cells
    i
    SAMPLES = cont_ca_sampler(Rsub(i,:),params);
    Rdcv(i,:) = make_mean_sample(SAMPLES,Rsub(i,:));
    S(i,:) = mean(samples_cell2mat(SAMPLES.ss,T));
end

temporal_time = toc(temp_time)/60

%%% Save results
if save_results
    cd(data_path)
    save(savename_CNMF_results,'A','R','Rsub','Rdcv','S','rows','cols')
end

total_time = toc(total_time)/60



