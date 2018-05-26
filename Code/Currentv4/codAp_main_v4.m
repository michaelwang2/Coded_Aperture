%% generate blur kernels
% codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% larger the window and blur kernel the better!
clear; close; clc;

load('codAp_v4_cropped_blur_patterns.mat');

nums = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
bdims = [33, 27, 33, 33, 37, 49, 49, 51, 57, 61];
KERS = cell(1, 10);
times = zeros(1, 10); 
dim = 500; % dimension of sharp and blurred images

for i = 5
    num = nums(i);
    bdim = bdims(i);
    [x, y] = size(Images{num});
    midx = round(x/2); midy = round(y/2);

    B = Images{num}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    S = Images{1}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    szKer = [bdim, bdim];
%     szFilt = [11, 11];
%     l1 = 10;
%     l2 = 3; % bigger less thick
%     l3 = l2;
%     l4 = 3; % bigger makes thicker
%     l5 = 0.0009; 
    tic
    % mKer = calcKer_quadprog(B, S, szKer, szFilt, l1, l2, l3, l4, l5);
    % nlevel = NoiseLevel(B);
    % l1/l2 = 40
    % l1 = 1/(nlevel^2);
    % l2 = 2*(2*floor(bdim/2) + 1)^2;
    l1 = 20;
    l2 = 8;
    mKer = calcKer_lsqnonneg(B, S, szKer, l1, l2);
    times(i) = toc;
    KERS{i} = mKer;
    fprintf('Kernel number %d is done!\n', i);
end

%% deconvolve, local energy estimate, and depth Map #1 v4 scene 2 sliding window flipped deconvL2 no weights
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
I = imread('codAp_v4_scene2.CR2');
I_cropped = I(300:850, 350:1150, :);
I_gray_double = im2double(rgb2gray(I_cropped));
[ry, cy] = size(I_gray_double);

% deconvL2 parameters
we = 0.006;
max_it = 200;

range = 2:10;
depthMapWindow = [25, 25];
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        % re-initialize
        minLocalEnergy = inf;
        minInd = 0;
        
        % obtain local window of scene
        I_local = I_gray_double((r-wingr):(r+wingr), (c-wingc):(c+wingc));
        for i = range
            % rotate kernel to match convention of deconSps %%%
            % tempDeconv = deconvSps(I_local, rot90(mKers{i}, 2), we, max_it);
            tempDeconv = deconvL2(I_local, rot90(KERS{i}, 2), we, max_it); % experient with rotation % Use the right orientation Kernel!!!!!!
            % tempDeconv = deconvL2(I_local, rot90(mKers{i}, 2), we);
            
            % calculate reconstruction error for local window
            if size(KERS{i}, 1) > size(tempDeconv, 1)
                reconv = conv2(KERS{i}, tempDeconv, 'valid');
            else
                reconv = conv2(tempDeconv, KERS{i}, 'valid');
            end
            reconv_wing = floor(size(reconv, 1)/2);
            reconError = I_local((wingr+1-reconv_wing):(wingr+1+reconv_wing), (wingc+1-reconv_wing):(wingc+1+reconv_wing)) - reconv;
            
            % calculate average local energy for this window (Frobenius
            % norm)
            % avgLocalEnergy = sum(sum(reconError.^2));
            avgLocalEnergy = mean(mean(reconError))^2;
            
            % find minimum energy
            if avgLocalEnergy < minLocalEnergy
                minLocalEnergy = avgLocalEnergy; % update
                minInd = i;
            end
        end
        % update depthMap 
        depthMap(r, c) = depths(minInd);
        disp(depths(minInd));
    end
    fprintf('Row: %i\n', r);
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

figure;
imshow(I_cropped);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% deconvolve test scene with all blur kernels #2 v4 scene 1
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
I = imread('codAp_v4_scene1.CR2');
I_cropped = I(300:850, 350:1150, :);

we = 0.006;
max_it = 200;

numBlurred = length(KERS);
deconvImages = cell(numBlurred, 1);
times = zeros(numBlurred, 1);
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
    deconvImages{i} = deconvL2(im2double(rgb2gray(I_cropped)), rot90(KERS{i}, 2), we, max_it); % Use the right orientation Kernel!!!!!!
    times(i) = toc;
    fprintf('Image %d\n', i);
end


%% local Energy Estimate and depth map #2 v4 scene 1 (weights only)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene1.mat');
load('w4_v3.mat');

numBlurred = length(KERS);
depthMapWindow = window;
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene1.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
range = 3:10;
weights = minw;
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
    localEnergyEst{i} = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
    for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
        for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
            avgLocalEnergy = sum(sum(reconErrors{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)).^2));
            localEnergyEst{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = avgLocalEnergy*weights(i);
        end
    end
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
depthMap = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
    for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
        minEng = inf;
        le_plot = zeros(1, length(range));
        for i = range
            le_plot(i-range(1)+1) = localEnergyEst{i}(r, c);
            curEng = localEnergyEst{i}(r, c);
            if curEng < minEng
                minEng = curEng;
                depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Local Energy','Rotation',0);
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% deconvolve test scene with all blur kernels #2 v4 scene 2
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
I = imread('codAp_v4_scene2.CR2');

% for deconvL2
we = 0.006;
max_it = 200;

% for deconvSps


numBlurred = length(KERS);
deconvImages = cell(numBlurred, 1);
times = zeros(numBlurred, 1);
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
    tdeconvImage = deconvL2(im2double(rgb2gray(I)), KERS{i}, we, max_it); % Use the right orientation Kernel!!!!!!
    deconvImages{i} = tdeconvImage(300:850, 350:1150);
    times(i) = toc;
    fprintf('Image %d\n', i);
end


%% local Energy Estimate and depth map #2 v4 scene 2 (weights only)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2.mat');
load('w4_v3.mat');

numBlurred = length(KERS);
depthMapWindow = window;
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% SSIM parameters
K1 = 0.000001;
K2 = 0.000001; 

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
% range = 3:10;
weights = minw;
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
    localEnergyEst{i} = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
    SSIMest{i} = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
    for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
        for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
            avgLocalEnergy = sum(sum(reconErrors{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)).^2));
            localEnergyEst{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = avgLocalEnergy*weights(i);
            SSIMest{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = ...
                SSIM(y(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)), conv2(deconvImages{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)), KERS{i}, 'same'), K1, K2);
        end
    end
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
depthMap = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
depthMapSSIM = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
    for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
        minEng = inf;
        maxSSIM = -inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            le_plot(i-range(1)+1) = localEnergyEst{i}(r, c);
            curEng = localEnergyEst{i}(r, c);
            if curEng < minEng
                minEng = curEng;
                depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
            
            curSSIM = SSIMest{i}(r, c);
            if curSSIM > maxSSIM
                maxSSIM = curSSIM;
                depthMapSSIM(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (weights only)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2_2.mat');
load('w7_25x25_deconvL2_2-10_nobias.mat');

numBlurred = length(KERS);
depthMapWindow = window; 
reconErrors = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
% range = 2:10;
weights = minw; 
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        minEng = inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            temp = reconErrors{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            curEng = sum(sum(temp.^2))*weights(i);
            le_plot(i-range(1)+1) = curEng;
            if curEng < minEng
                minEng = curEng;
                depthMap(r, c) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (weights only) (flipped kernel for deconvL2)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2_2_flipped.mat');
load('w8_19x19_deconvL2_2-10_nobias_flipped.mat');

numBlurred = length(KERS);
depthMapWindow = window; 
reconErrors = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
weights = minw; 
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        minEng = inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            temp = reconErrors{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            curEng = sum(sum(temp.^2))*weights(i);
            le_plot(i-range(1)+1) = curEng;
            if curEng < minEng
                minEng = curEng;
                depthMap(r, c) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 (bias + weights)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2.mat');
load('w4_v5.mat');

numBlurred = length(KERS);
depthMapWindow = window; 
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% SSIM parameters
K1 = 0.000001;
K2 = 0.000001;

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
% range = 2:10;
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
    localEnergyEst{i} = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
    SSIMest{i} = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
    for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
        for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
            avgLocalEnergy = sum(sum(reconErrors{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)).^2));
            localEnergyEst{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = avgLocalEnergy*minw(i) + minb(i);
            SSIMest{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = ...
                SSIM(y(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)), conv2(deconvImages{i}(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)), KERS{i}, 'same'), K1, K2);
        end
    end
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
depthMap = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
depthMapSSIM = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
    for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
        minEng = inf;
        maxSSIM = -inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            le_plot(i-range(1)+1) = localEnergyEst{i}(r, c);
            curEng = localEnergyEst{i}(r, c);
            if curEng < minEng
                minEng = curEng;
                depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
            
            curSSIM = SSIMest{i}(r, c);
            if curSSIM > maxSSIM
                maxSSIM = curSSIM;
                depthMapSSIM(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;
% figure;
% imagesc(depthMapSSIM);
% axis equal;
% title('Depth Map');
% colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (weights and bias)
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2_2_flipped.mat');
load('w8_19x19_deconvL2_2-10_bias_flipped.mat');

numBlurred = length(KERS);
depthMapWindow = window;  
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

% range = 2:numBlurred; % range of blur kernels using
% weights = linspace(1, 1, numBlurred);
% range = 2:9;
weights = minw; 
bias = minb;
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        minEng = inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            temp = reconErrors{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            curEng = sum(sum(temp.^2))*weights(i) + bias(i);
            le_plot(i-range(1)+1) = curEng;
            if curEng < minEng
                minEng = curEng;
                depthMap(r, c) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (no weights or biases) deconvL2
clear; close; clc;
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvL2_scene2_2_flipped.mat');

numBlurred = length(KERS);
depthMapWindow = [31, 31];  
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

range = 2:10;
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, KERS{i}, 'same');
end

% generate depth map
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        minEng = inf;
        le_plot = zeros(1, length(range)); 
        for i = range
            temp = reconErrors{i}((r-wingr):(r+wingr), (c-wingc):(c+wingc));
            curEng = sum(sum(temp.^2));
            le_plot(i-range(1)+1) = curEng;
            if curEng < minEng
                minEng = curEng;
                depthMap(r, c) = depths(i);
            end
        end
        % plot(depths(range), le_plot);
    end
end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));

grid on; box on;
xlabel('Depth [m]');
yy = ylabel('Error');
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

figure;
imshow(y_orig);
figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;

%% local Energy Estimate and depth map #2 v4 scene 2 sliding window (no weights or biases) deconvSps
clear; close; clc;
addpath('Kernels');
load('Kernels/codAp_v4_bk_cropped.mat');
load('codAp_v4_deconvSps_scene2_flipped_rgb.mat');
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
depthMapWindow = [25, 25]; 
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y_orig = imread('codAp_v4_scene2.CR2');
y_orig = y_orig(300:850, 350:1150, :);
y = im2double(rgb2gray(y_orig(:,:,:)));
[ry, cy] = size(y);

range = 3:10;
y_other = cell(numBlurred, 1);
for i = range
    y_other{i} = conv2(deconvImages{i}, KERS{i}, 'same');
    reconErrors{i} = y - y_other{i};
end

% generate depth map
K1 = 1e-3; K2 = 1e-3;
figure; hold on;
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m
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


