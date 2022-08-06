%%% Overlay the good ROIs onto static ims to classify cells as rabies-positive or not
close all
clear all

%parameters
PC_choice_new = 3; % corresponds to 7 PCs
PC_choice_old = 2; % corresponds to 6 PCs
CorS = 2; %Calcium or Spiking: 1 = Ca, 2 = Sp
orig_line_leng = 256; %In pixels
grn_stack_prefix = '_Cycle00001_Ch3_0000';
red_stack_prefix = '_Cycle00001_Ch2_0000';
pix_val_lb = 0.005; %quantile of pixel vals below which to ignore
pix_val_ub = 0.995; %quantile of pixel vals above which to set to 1
us_f = 100;

%Initialize
next_fig_num = 1;
roi_max_thr = 0.25;
sim_pre_post = zeros(9,2);

%Define files to use
filename(1).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(2).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(3).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(4).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(5).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(6).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D08';
filename(7).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D16';
filename(8).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D16';
filename(9).ts_dir = 'C:\Users\jam5064\Desktop\Current Data\D16';

filename(1).align_old = 'D08_10272017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(2).align_old = 'D08_10292017_align_v12_001_CNMF_results_1nb_1_5_p1.mat';
filename(3).align_old = 'D08_10302017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(4).align_old = 'D08_11012017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(5).align_old = 'D08_11022017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(6).align_old = 'D08_11032017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(7).align_old = 'D16_11162017_align_v12_001_CNMF_results_1nb_1_5_p1.mat';
filename(8).align_old = 'D16_11172017_align_v12_001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(9).align_old = 'D16_11222017b_align_v12_001_CNMF_results_1nb_1_5_p1.mat';


filename(1).align_new = 'D08_10272017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(2).align_new = 'D08_10292017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1.mat';
filename(3).align_new = 'D08_10302017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(4).align_new = 'D08_11012017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(5).align_new = 'D08_11022017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(6).align_new = 'D08_11032017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(7).align_new = 'D16_11162017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1.mat';
filename(8).align_new = 'D16_11172017_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(9).align_new = 'D16_11222017b_align_AIQ_75_off_v3001_CNMF_results_1nb_1_5_p1.mat';


filename(1).mi = 'D08_10272017_mean_ims_128.mat';
filename(2).mi = 'D08_10292017_mean_ims_128.mat';
filename(3).mi = 'D08_10302017_mean_ims_128.mat';
filename(4).mi = 'D08_11012017_mean_ims_128.mat';
filename(5).mi = 'D08_11022017_mean_ims_128.mat';
filename(6).mi = 'D08_11032017_mean_ims_128.mat';
filename(7).mi = 'D16_11162017_mean_ims_128.mat';
filename(8).mi = 'D16_11172017_mean_ims_128.mat';
filename(9).mi = '2205_D16_11222017_mean_ims_128.mat';


filename(1).cnmf = 'TSeries-10272017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(2).cnmf = 'TSeries-10292017-D08-001_CNMF_results_1nb_1_5_p1.mat';
filename(3).cnmf = 'TSeries-10302017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(4).cnmf = 'TSeries-11012017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(5).cnmf = 'TSeries-11022017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(6).cnmf = 'TSeries-11032017-D08-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(7).cnmf = 'TSeries-11162017-D16-001_CNMF_results_1nb_1_5_p1.mat';
filename(8).cnmf = 'TSeries-11172017-D16-001_CNMF_results_1nb_1_5_p1_fin.mat';
filename(9).cnmf = 'TSeries-11222017b-D16-001_CNMF_results_1nb_1_5_p1.mat';


filename(1).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-005';
filename(1).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-006';

filename(2).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-003';
filename(2).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-004';

filename(3).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-008';
filename(3).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-009';

filename(4).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-010';
filename(4).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-011';

filename(5).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-012';
filename(5).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-013';

filename(6).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-014';
filename(6).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D08\ZSeries\ZSeries-11112017-D08-015';

filename(7).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-001';
filename(7).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-002';

filename(8).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-003';
filename(8).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-004';

filename(9).grn_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-011';
filename(9).red_dir = 'C:\Users\jam5064\Desktop\Current Data\D16\ZSeries\ZSeries-11302017-D16-012';


filename(1).z_length = 41;
filename(2).z_length = 41;
filename(3).z_length = 41;
filename(4).z_length = 41;
filename(5).z_length = 41;
filename(6).z_length = 41;
filename(7).z_length = 41;
filename(8).z_length = 41;
filename(9).z_length = 41;


filename(1).clsfd = 'D08_10272017_classified.mat';
filename(2).clsfd = 'D08_10292017_classified.mat';
filename(3).clsfd = 'D08_10302017_classified.mat';
filename(4).clsfd = 'D08_11012017_classified.mat';
filename(5).clsfd = 'D08_11022017_classified.mat';
filename(6).clsfd = 'D08_11032017_classified.mat';
filename(7).clsfd = 'D16_11162017_classified.mat';
filename(8).clsfd = 'D16_11172017_classified.mat';
filename(9).clsfd = 'D16_11222017b_classified.mat';


%Main FOR loop
run_time = tic;
for i = 1:9
    i
    
    %Load data
    cd(filename(i).ts_dir)
        
    load(filename(i).cnmf);
    load(filename(i).mi);
    if ~isempty(filename(i).clsfd)
        load(filename(i).clsfd)
    end
    load(filename(i).align_new,'good_cells_full')
    if CorS == 1
        good_cells_new = good_cells_full(PC_choice_new).ca;
    else good_cells_new = good_cells_full(PC_choice_new).sp;
    end
    num_cells = length(good_cells_new)
 
     
    row_start =  sum(rows < orig_line_leng/2) + 1;
    col_start = sum(cols < orig_line_leng/2) + 1;
    
    mean_mov(1).im = mean_im_start;
    mean_mov(2).im = mean_im_mid;
    mean_mov(3).im = mean_im_end;
    mms = size(mean_mov(1).im);

    %Now collect with ROIs into an array - both scaled and binary
    num_rows = orig_line_leng - length(rows);
    num_cols = orig_line_leng - length(cols);

    ROIs = zeros(num_rows,num_cols,num_cells);
    for k = 1:num_cells
        current_cell = good_cells_new(k);
        ROIs(:,:,k) = reshape(A(:,good_cells_new(k)),num_rows,num_cols);
    end

    %Make full size version of ROIs for potential ease later
    ROIs_full = zeros(orig_line_leng,orig_line_leng,num_cells);
    for k = 1:num_cells
        ROIs_full(row_start:row_start+mms(1)-1,col_start:col_start+mms(2)-1,k) =  ROIs(:,:,k);
    end
    ROIs_full_sum = sum(ROIs_full,3);
    A_full = reshape(ROIs_full,orig_line_leng*orig_line_leng,num_cells);
    A_full_sparse = sparse(A_full);
    
    %Find best overlay with mean_im given ROIs
    sim_overlays = zeros(3,1);
    for k = 1:3
        mean_mov_full = zeros(orig_line_leng,orig_line_leng);
        mean_mov_full(row_start:row_start+mms(1)-1,col_start:col_start+mms(2)-1) = mean_mov(k).im;
        sim_overlays(k) = corr2(ROIs_full_sum,mean_mov_full);
    end
    [~,ind] = max(sim_overlays)

    mean_mov = mean_mov(ind).im;
    mean_mov_full = zeros(orig_line_leng,orig_line_leng);
    mms = size(mean_mov);
    mean_mov_full(row_start:row_start+mms(1)-1,col_start:col_start+mms(2)-1) = mean_mov;
       
    %sanity check - compare ROI map to mean_mov
    figure(6)
    imagesc(mean_mov,[min(mean_mov(:)),max(mean_mov(:))]);
    axis tight; axis equal;
    hold on;
    for k = 1:num_cells
        roi_temp = medfilt2(ROIs(:,:,k),[2,2]);
        roi_temp(roi_temp<roi_max_thr*max(roi_temp(:))) = 0;
        BW = bwareafilt(roi_temp>0,1);                
        BW2 = bwboundaries(BW);
        cent = regionprops(BW,'centroid');
        for ii = 1:length(BW2)
            BW2{ii} = fliplr(BW2{ii});
            plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color','m', 'linewidth', 2);
        end
        text(cent.Centroid(1)-8,cent.Centroid(2)-8,strtrim(cellstr(num2str(k))),'color',[1,.5 ,0],'fontsize',16,'fontname','helvetica','fontweight','bold');
    end
    title('ROI overlay - looks right?')
    set(gcf,'Position',[300 0 600 500])
     
    % Determine which green stack image to use with a nonrigid correction 
    % Compare two methods now
    % Find the maximal xcorr between the mean_mov and the different elements of the stack
    sim_vals = zeros(filename(i).z_length,1); 
    sharpness_array = zeros(filename(i).z_length,2);
    for k = 1:filename(i).z_length
        k
        cd(filename(i).grn_dir)
        names = strsplit(filename(i).grn_dir,'\');
        grn_stk_im = loadtiff([char(names{end}) grn_stack_prefix num2str(k,'%02.f') '.ome.tif']);

        [~,grn_stk_im] = correct_bidirectional_offset(grn_stk_im,1,10);
        lb_gs = double(quantile(grn_stk_im(:),pix_val_lb));
        ub_gs = double(quantile(grn_stk_im(:),pix_val_ub));
        gs_g = mat2gray(grn_stk_im, [lb_gs ub_gs]);
        
        %Method One: NormCorre
        corr_time_1 = tic;
        options_nonrigid = NoRMCorreSetParms('d1',size(mean_mov_full,1),'d2',size(mean_mov_full,2),'grid_size',[32,32],'overlap_pre',[16,16],'max_shift',30,'max_dev',[20,20],'us_fac',us_f,'upd_template',0);
        [mean_mov_corr_1,~,~] = normcorre(mean_mov_full,options_nonrigid,grn_stk_im);
        t1 = toc(corr_time_1)
        sim_vals(k) = corr2(mean_mov_corr_1,grn_stk_im);
        
        % Estimate sharpness using the gradient magnitude.
        % sum of all gradient norms / number of pixels give us the sharpness metric.
        [Gx, Gy] = gradient(gs_g);
        S = sqrt(Gx.*Gx+Gy.*Gy);
        sharpness_array(k,1) = sum(sum(S))./(numel(Gx));
        
        cd(filename(i).red_dir)
        names = strsplit(filename(i).red_dir,'\');
        red_stk_im = loadtiff([char(names{end}) red_stack_prefix num2str(k,'%02.f') '.ome.tif']);

        [~,red_stk_im] = correct_bidirectional_offset(red_stk_im,1,10);
        lb_rs = double(quantile(red_stk_im(:),pix_val_lb));
        ub_rs = double(quantile(red_stk_im(:),pix_val_ub));
        rs_g = mat2gray(red_stk_im, [lb_rs ub_rs]);

        [Gx, Gy] = gradient(rs_g);
        S = sqrt(Gx.*Gx+Gy.*Gy);
        sharpness_array(k,2) = sum(sum(S))./(numel(Gx));
    end
    
    % Sanity check - look at the similarity vals to make sure they look sensible
    figure(10)
    plot(1:filename(i).z_length,sim_vals)
    xlabel('Stack number')
    ylabel('image similarity - Method 1')
    
    % Find max overlap; make sure images aren't particularly blurry
    [sorted,inds] = sort(sim_vals,'descend');
    best_stk_im = inds(1);
    
%     if i == 18
%         best_stk_im = 25;
%     end
    
%     %Check if the best green and red are blurry from movement
%     current_rank = 1;
%     flag = 0;
%     while flag == 0
%         if sharpness_array(inds(current_rank),1) < mean(sharpness_array(:,1))-3*std(sharpness_array(:,1)) || sharpness_array(inds(current_rank),2) < mean(sharpness_array(:,2))-3*std(sharpness_array(:,2))
%             i
%             disp('image was blurry, trying next best')
%             current_rank = current_rank + 1
%         else
%             flag = 1;
%         end
%     end
%     best_stk_im = inds(current_rank)
%     
%     if current_rank > 1
%         pause %pause to assess if everything is copacetic
%     end
       
    % Load the green and red stack images to be used for the rest of the analysis
    cd(filename(i).grn_dir)
    names = strsplit(filename(i).grn_dir,'\');
    grn_stk_im = loadtiff([char(names{end}) grn_stack_prefix num2str(best_stk_im,'%02.f') '.ome.tif']);
    
    [~,grn_stk_im] = correct_bidirectional_offset(grn_stk_im,1,10);
    lb_gs = double(quantile(grn_stk_im(:),pix_val_lb));
    ub_gs = double(quantile(grn_stk_im(:),pix_val_ub));
    gs_g = mat2gray(grn_stk_im, [lb_gs ub_gs]);
    [gs_ind, ~] = gray2ind(gs_g,256);
    filler = zeros(size(gs_ind),'uint8');
    gsImage = cat(3,filler,gs_ind,filler);
    figure(1)
    imshow(gsImage);
    title('Green Static Image (grn)')
    axis tight; axis equal
    set(gcf,'Position',[0 0 600 500])
    
    cd(filename(i).red_dir)
    names = strsplit(filename(i).red_dir,'\');
    red_stk_im = loadtiff([char(names{end}) red_stack_prefix num2str(best_stk_im,'%02.f') '.ome.tif']);
    
    [~,red_stk_im] = correct_bidirectional_offset(red_stk_im,1,10);
    lb_rs = double(quantile(red_stk_im(:),pix_val_lb));
    ub_rs = double(quantile(red_stk_im(:),pix_val_ub));
    rs_g = mat2gray(red_stk_im, [lb_rs ub_rs]);
    [rs_ind, ~] = gray2ind(rs_g,256);

%     %%Sanity check: plot red and green overlay pre and post shift to compare
%     rgImage = cat(3,rs_ind,gs_ind,filler);
%     figure(next_fig_num); next_fig_num = next_fig_num+1;
%     imshow(rgImage);
%     title('Green Stack Image (grn) vs. Red Stack Image (red)')
%     axis image
    
    %[yoffset,xoffset,max_corr]=findoff(orig,shifted) - positive values mean shifted is offset down and to the right
    [yoffset,xoffset,~] = findoff(gs_g,rs_g)

    % Trim or expand to sensibly deal with the offsets (verified this chunk works with sample matrices, pos and neg shifts, when shifts are same magnitude
    red_stk_im_shifted = red_stk_im;
    if yoffset > 0
        red_stk_im_shifted = [red_stk_im_shifted(yoffset+1:end,:); zeros(yoffset,orig_line_leng)];
    elseif yoffset < 0
        red_stk_im_shifted = [zeros(abs(yoffset),orig_line_leng); red_stk_im_shifted(1:end-abs(yoffset),:)];
    end

    if xoffset > 0
        red_stk_im_shifted = [red_stk_im_shifted(:,xoffset+1:end) zeros(orig_line_leng,xoffset)];
    elseif xoffset< 0
        red_stk_im_shifted = [zeros(orig_line_leng,abs(xoffset)) red_stk_im_shifted(:,1:end-abs(xoffset))];
    end

    %Plot the shifted overlay
    lb_red = double(quantile(red_stk_im_shifted(:),pix_val_lb));
    ub_red = double(quantile(red_stk_im_shifted(:),pix_val_ub));
    rss_g = mat2gray(red_stk_im_shifted, [lb_red ub_red]); 
    [rss_ind, ~] = gray2ind(rss_g,256);

    rgImage_gs_rss = cat(3,rss_ind,gs_ind,filler);
    figure(2)
    imshow(rgImage_gs_rss);
    title('Green Stack Image (grn) vs. SHIFTED Red Stack Image (red)')
    axis tight; axis equal
    set(gcf,'Position',[0 100 600 500])
    
    rImage = cat(3,rss_ind,filler,filler);
    figure(3)
    imshow(rImage);
    title('SHIFTED Red Stack Image (red)')
    axis tight; axis equal
    set(gcf,'Position',[200 0 600 500])
    
    % Normcorre-based alignment
    im_sim_pre = corr2(mean_mov_full,grn_stk_im);
    corr_time_1 = tic;
    options_nonrigid = NoRMCorreSetParms('d1',size(mean_mov_full,1),'d2',size(mean_mov_full,2),'grid_size',[32,32],'overlap_pre',[16,16],'max_shift',30,'max_dev',[20,20],'us_fac',us_f,'upd_template',0);
    [mean_mov_shifted,align_shifts,~] = normcorre(mean_mov_full,options_nonrigid,grn_stk_im);
    toc(corr_time_1)
    im_sim_post = corr2(mean_mov_shifted,grn_stk_im);
    sim_pre_post(i,:) = [im_sim_pre im_sim_post];
    
    lb_mms = double(quantile(mean_mov_shifted(:),pix_val_lb)); %dont include zero padding into lb and ub calculation
    ub_mms = double(quantile(mean_mov_shifted(:),pix_val_ub));
    mms_g = mat2gray(mean_mov_shifted, [lb_mms ub_mms]);
    [mms_ind, ~] = gray2ind(mms_g,256);
    mmsImage = cat(3,filler,mms_ind,filler);
%     figure(next_fig_num); next_fig_num = next_fig_num+1;
%     imshow(mmsImage);
%     title('SHIFTED Mean Movie Image (grn)')
%     axis image
    
    %Alignment method #1: Use shift reconstruct to map the ROIs to proper space
%     I = shift_reconstruct(Y,shifts,diffphase,us_fac,Nr,Nc,Np,method,add_value)
    shifts_fov = reshape(imresize(align_shifts(1).shifts,[orig_line_leng,orig_line_leng]),[],2);
    shifts_components = sparse(diag(1./sum(A_full_sparse)))*A_full_sparse'*shifts_fov;
    A_full_shifted = 0*A_full_sparse;
    for j = 1:num_cells
       a_temp = reshape(full(A_full_sparse(:,j)),[orig_line_leng,orig_line_leng,1]);
       a_temp = shift_reconstruct(a_temp,shifts_components(j,:),0);   
       A_full_shifted(:,j) = sparse(a_temp(:));
    end
    
    ROIs_full_shifted = zeros(orig_line_leng,orig_line_leng,num_cells);
    for k = 1:num_cells
        ROIs_full_shifted(:,:,k) = reshape(A_full_shifted(:,k),orig_line_leng,orig_line_leng);
    end
    
    %sanity check - compare ROI map to mean_mov
%     figure(next_fig_num); next_fig_num = next_fig_num+1;
    figure(4)
    imagesc(mean_mov_shifted,[min(mean_mov(:)),max(mean_mov(:))]); %dont include zero padding into lb and ub calculation
    axis tight; axis equal;
    hold on;
    centroids = zeros(num_cells,2);
    for k = 1:num_cells
        roi_temp = medfilt2(ROIs_full_shifted(:,:,k),[2,2]);
        roi_temp(roi_temp<roi_max_thr*max(roi_temp(:))) = 0;
        BW = bwareafilt(roi_temp>0,1);  
        cent = regionprops(BW,'centroid');
        BW2 = bwboundaries(BW);
        for ii = 1:length(BW2)
            BW2{ii} = fliplr(BW2{ii});
            plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color','m', 'linewidth', 2);
        end
        text(cent.Centroid(1)-8,cent.Centroid(2)-8,strtrim(cellstr(num2str(k))),'color',[1,.5 ,0],'fontsize',16,'fontname','helvetica','fontweight','bold');
        centroids(k,:) = [cent.Centroid(1) cent.Centroid(2)];
    end
    title('SHIFTED ROI overlay - looks right?')
    axis tight; axis equal;
    set(gcf,'Position',[100 100 600 500])

    figure(5)
    imshow(rgImage_gs_rss);
    title('SHIFTED ROI overlay on Green Stack Image (grn) vs. SHIFTED Red Stack Image (red)')
    axis tight; axis equal;
    hold on;
    for k = 1:num_cells
        roi_temp = medfilt2(ROIs_full_shifted(:,:,k),[2,2]);
        roi_temp(roi_temp<roi_max_thr*max(roi_temp(:))) = 0;
        BW = bwareafilt(roi_temp>0,1);                
        BW2 = bwboundaries(BW);
        for ii = 1:length(BW2)
            BW2{ii} = fliplr(BW2{ii});
            plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color','m', 'linewidth', 2);
        end
    end
    axis tight; axis equal;
    set(gcf,'Position',[200 100 600 500])
    
    cells2classify = cells2classify  %to show me what's needed
    for k = 1:length(cells2classify)
        kk = cells2classify(k)
        centroid = centroids(kk,:)
        xx = input('3 = def labeled, 2 = probably labeled, 1 = def unlabeled, 0 = unclear') ; 
        if isempty(xx)
            classification_new(kk,2) = 0;
            disp('hit enter, returned "undefined" ')
        else
            classification_new(kk,2) = xx;
        end
    end
    cd(filename(i).ts_dir)
    names = strsplit(filename(i).cnmf,'-');
    names2 = strsplit(names{4},'_');
    savename = [names{3} '_' names{2} '_' names2{1} '_classified_v3'];
    classification = classification_new;
    save(savename,'classification')
    disp([savename ' saved'])
    
    classification
    disp(['Make sure it looks right'])
    pause    
    
    clear mean_mov
    clear classification_new
    clear classification
    clear cells2classify
end
toc(run_time)       



