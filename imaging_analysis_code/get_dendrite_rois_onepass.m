function [binaryImageNew, centroids] = get_dendrite_rois_onepass(meanMovIm)

pix_val_lb = 0.005; %quantile of pixel vals below which to ignore
pix_val_ub = 0.995; %quantile of pixel vals above which to set to 1
captionFontSize = 14;

%find bounds and plot
lb_mm = double(quantile(meanMovIm(:),pix_val_lb));
ub_mm = double(quantile(meanMovIm(:),pix_val_ub));
figure; imagesc(meanMovIm,[lb_mm ub_mm]); axis image

%convert to ind
originalImage = mat2gray(meanMovIm, [lb_mm ub_mm]);
[originalImage, ~] = gray2ind(originalImage,256);

% histogram of intensities
[pixelCount, grayLevels] = imhist(originalImage);
figure;
bar(pixelCount);
title('Histogram of original image', 'FontSize', captionFontSize);

%start to segment
thresholdValue = 110;
binaryImage = originalImage > thresholdValue; % Bright objects will be chosen if you use >.

% Do a "hole fill" to get rid of any background pixels or "holes" inside the blobs.
binaryImage = imfill(binaryImage, 'holes');
figure;
imshow(binaryImage); 
title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize); 

%Now try the watershed on the distance transform on the thresholded images
D = bwdist(~binaryImage);
D = -D;
D(~binaryImage) = Inf;
L = watershed(D);
L(~binaryImage) = 0;
rgb = label2rgb(L,'jet',[.5 .5 .5]);
figure
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of original with distance transform')

%get rid of elements that are too small - less than 4 pixels
num_cells = max(L(:));
for i = 1:num_cells
    num_pixels = length(find(L(:)==i));
    if num_pixels < 8
        L(L==i) = 0;
    end
end
binaryImageNew = L > 0.5;
L = bwlabel(binaryImageNew);
L(~binaryImageNew) = 0;
num_cells = max(L(:));
rgb = label2rgb(L,'jet',[.5 .5 .5]);
figure
imshow(rgb,'InitialMagnification','fit')
title(['After removing small cells, num cell = ' num2str(num_cells)])

% % bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% % Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.
% figure;
% imshow(originalImage);
% title('Outlines, from bwboundaries()', 'FontSize', captionFontSize); 
% axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
% hold on;
% boundaries = bwboundaries(binaryImageNew);
% numberOfBoundaries = size(boundaries, 1);
% for k = 1 : numberOfBoundaries
% 	thisBoundary = boundaries{k};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 0.5);
% end
% set(gcf,'Position',[10 10 1600 1200]);
% hold off;

stats = regionprops(bwlabel(binaryImageNew),'centroid');
centroids = cat(1, stats.Centroid);

    

