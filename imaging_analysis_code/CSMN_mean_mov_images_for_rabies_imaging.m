close all
clear all

filenames(1).mat = 'TSeries-11012017-D19-001_mov_ch3_mov_nonrigid_corr.mat';
% filenames(2).mat = 'TSeries-11022017-D20-002_mov_ch3_mov_nonrigid_corr.mat';
% filenames(3).mat = 'TSeries-11032017-D20-001_mov_ch3_mov_nonrigid_corr.mat';
% filenames(4).mat = 'TSeries-11042017-D20-001_mov_ch3_mov_nonrigid_corr.mat';
% filenames(5).mat = 'TSeries-11052017-D20-001_mov_ch3_mov_nonrigid_corr.mat';
% % filenames(6).mat = 'TSeries-11052017-D19-002_mov_ch3_mov_nonrigid_corr.mat';

frames2use = 100;
using_corrected = 1;
lb = 0.05; %quantile of pixel vals below which to ignore
ub = 0.995; %quantile of pixel vals above which to set to 1

for i = 1:length(filenames)
    m = matfile(filenames(i).mat);
    if using_corrected
        Y = m.mov_nonrigid_corr(:,:,end-frames2use:end);
    else Y = m.mov(:,:,end-frames2use:end);
        clear mov
    end
%     if using_corrected
%         Y = m.mov_nonrigid_corr(:,:,1:frames2use);
%     else Y = m.mov(:,:,1:frames2use);
%         clear mov
%     end
    names = strsplit(filenames(i).mat,'-');

    lb_gs = double(quantile(Y(:),lb));
    ub_gs = double(quantile(Y(:),ub));
    gr = mat2gray(mean(Y,3), [lb_gs ub_gs]); 
    [gr_ind, cm] = gray2ind(gr,256);
    filler = zeros(size(gr_ind),'uint8');
    grnImage = cat(3,filler,gr_ind,filler);
    figure;
    imshow(grnImage);
    title([names(2) ' ' names(3)])
    
    B = imgaussfilt(grnImage,0.75);
    figure;
    imshow(B);

end
        
      

