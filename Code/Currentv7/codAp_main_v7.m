%% generate blur kernels
% codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% larger the window and blur kernel the better!
clear; close; clc;

load('codAp_v7_1_4_cropped_resized_blur_patterns.mat');

nums = [2, 3, 4, 5, 6, 7, 8, 9, 10];
bdims = [33, 27, 33, 33, 37, 49, 49, 51, 57];
KERS = cell(1, 9);
times = zeros(1, 9); 
dim = 500; % dimension of sharp and blurred images

for i = 1:9
    num = nums(i);
    bdim = bdims(i);
    [x, y] = size(Images{num});
    midx = round(x/2); midy = round(y/2);

    B = Images{num}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    S = Images{1}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    szKer = [bdim, bdim];
    tic;
    szFilt = [11, 11];
    l1 = 10;
    l2 = 0; % bigger less thick
    l3 = 0;
    l4 = 20; % bigger makes thicker
    l5 = 0; 
    mKer = calcKer_quadprog(B, S, szKer, szFilt, l1, l2, l3, l4, l5);
    % nlevel = NoiseLevel(B);
    % l1/l2 = 40
    % l1 = 1/(nlevel^2);
    % l2 = 2*(2*floor(bdim/2) + 1)^2;
%     l1 = 20;
%     l2 = 8;
%     mKer = calcKer_lsqnonneg(B, S, szKer, l1, l2);
    times(i) = toc;
    KERS{i} = mKer;
    fprintf('Kernel number %d is done!\n', i);
end

%% deconvolve test scene with all blur kernels #2 v4 scene 2 (color channels)
clear; close; clc;
% load('codAp_v7_bks_cropped.mat');
% I = imread('1_4_focus/scene1.CR2');
load('codAp_v6_1_3_bks.mat');
I = imread('v6_scene1.CR2');
% orig_crop = [500,1850,300,2900];
orig_crop = [550,1850,500,2700];

% for deconvL2
we = 0.006;
max_it = 200;

% for deconvSps


numBlurred = length(KERS);
deconvImages = cell(numBlurred, 3);
times = zeros(numBlurred, 1);
% MATLAB is rgb (1:3)
for c = 1:3
    fprintf('Color %d\n', c);
    for i = 1:numBlurred
        tic;
        % rotate kernel to match convention of deconSps %%%
        % deconvolve window by window so the artifacts do not affect other
        % scenes
        %{
        To obtain depth for the full image, we consider each sub-window independently, 
        trying a range of defocus scales and picking the one which minimizes the deconvolution error. 
        However, this can be noisy. 
        Also, like most passive depth estimation methods, our approach requires texture, 
        so has trouble with uniform areas of the image.
        Therefore, we use a markov random field to regularize the local depth map. 
        This prefers the depth to be piece-wise constant and that depth discontinuities should align with image discontinuities. 
        This gives the improved depth map on the right. 
        The colors correspond to depth from the camera, red being closer and blue being further away.
        %}
        % deconvImages{i} = deconvSps(im2double(rgb2gray(I_cropped)),rot90(mKers{i}, 2),we,max_it);
        tdeconvImage = deconvL2(im2double(I(:, :, c)), rot90(KERS{i}, 2), we, max_it); % Use the right orientation Kernel!!!!!!
        deconvImages{i, c} = tdeconvImage(orig_crop(1):orig_crop(2), orig_crop(3):orig_crop(4));
        times(i) = toc;
        fprintf('Image %d\n', i);
    end
end


%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (no weights or biases) deconvSps
clear; close; clc;
% load('codAp_v7_bks_cropped.mat');
% load('codAp_v7_deconvL2_scene1_flipped_rgb.mat');
load('codAp_v6_1_3_bks.mat');
load('codAp_v6_1_3_deconvL2_scene1_flipped.mat');
% load('codAp_v4_deconvL2_scene2_2_flipped.mat');
% load('w8_19x19_deconvL2_2-10_nobias_flipped.mat');

% color = 3; % rgb green > red
tempDeconv = deconvImages;
deconvImages = cell(10, 1);
for i = 1:10
    temprgb = zeros(size(tempDeconv{1,1}, 1), size(tempDeconv{1,1}, 2), 3);
    for c = 1:3
        temprgb(:, :, c) = tempDeconv{i, c};
    end
    deconvImages{i} = rgb2gray(temprgb);
end

numBlurred = length(KERS); 
depthMapWindow = [61, 61]; 
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
% y_orig = imread('1_4_focus/scene1.CR2');
y_orig = imread('v6_scene1.CR2');
y_orig = y_orig(orig_crop(1):orig_crop(2), orig_crop(3):orig_crop(4), :);
y = im2double(rgb2gray(y_orig(:,:,:)));
[ry, cy] = size(y);

range = 4:10;
y_other = cell(numBlurred, 1);
for i = range
    y_other{i} = conv2(deconvImages{i}, KERS{i}, 'same');
    reconErrors{i} = y - y_other{i};
end

% generate depth map
K1 = 1e-3; K2 = 1e-3;
figure; hold on;
% depths = 1.5:0.1:2.3; % m
depths = 1.4:0.1:2.3;
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
depthMapSSIM = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        minEng = inf;
        maxSSIM = -inf; 
        le_plot = zeros(1, length(range)); 
        for i = range
            y_patch = y((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            y_other_patch = y_other{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            
            temp = reconErrors{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            curEng = sum(sum(temp.^2));
            curSSIM = SSIM(y_patch, y_other_patch, K1, K2);
            le_plot(i-range(1)+1) = curEng;
            if curEng < minEng
                minEng = curEng;
                depthMap(r, c) = depths(i);
            end
            if curSSIM > maxSSIM
                maxSSIM = curSSIM;
                depthMapSSIM(r, c) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));
depthMapSSIM = depthMapSSIM((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map (Reconstruction Error)');
h = colorbar;
y = ylabel(h, 'Depth [m]', 'rotation', -90);
set(y, 'Units', 'Normalized', 'Position', [5, 0.5, 0]);

figure; 
imagesc(depthMapSSIM);
axis equal;
title('Depth Map (SSIM)');
h = colorbar;
y = ylabel(h, 'Depth [m]', 'rotation', -90);
set(y, 'Units', 'Normalized', 'Position', [5, 0.5, 0]);

% generate all focused image (grey)
all_focused_image = zeros(size(depthMap, 1), size(depthMap, 2), 3);
for i = 1:size(all_focused_image, 1)
    for j = 1:size(all_focused_image, 2)
        % figure out which depth 
        depth = depthMap(i, j);
        indd = find(depths == depth);
        for c = 1:3
            all_focused_image(i, j, c) = tempDeconv{indd, c}(wingr+i, wingc+j); % deconvImagesrgb{indd, c}(wingr+i, wingc+j);
        end
    end
end
figure;
imshow(all_focused_image);


