% load D
Y = mov_nonrigid_corr(:,:,1:1000);
m = min(Y(:));
Y = Y-m;
nnY = double(quantile(Y(:),0.005));
mmY = double(quantile(Y(:),0.995));

s = size(Y);
Y_ind = zeros(s(1),s(2),1,s(3));
for i = 1:s(3)
    gm = mat2gray(Y(:,:,i), [nnY mmY]); 
    [Y_ind(:,:,:,i), cm] = gray2ind(gm,256);
end
mov = immovie(Y_ind+1,cm); %%%  X is an m-by-n-by-1-by-k array, where k is the number of images.
implay(mov)


v = VideoWriter('D12_0928_001.avi','Uncompressed AVI');
open(v)
writeVideo(v,mov)
close(v)