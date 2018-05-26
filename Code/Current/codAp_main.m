clear; close; clc;
% Toggles
use_undistorted_images = true; % true or false
aperture = 'coded'; % conventional or coded
numBlurred = 10; % number of blurred calibration images
scene = 1; % scene number to use for depth generation and all-focus image
blurKerDims = [11, 11]; % dimension of the desired blur kernels
depthMapWindow = [9, 9]; % dimension of a window of constant depth 

% add paths
addpath(genpath('Camera Calibration'));
addpath(genpath('Depth Generation'));

% file prefix
if use_undistorted_images
    if strcmp(aperture, 'coded')
        filePrefix = 'Coded_Undist';
    elseif strcmp(aperture, 'conventional')
        filePrefix = 'Conventional_Undist';
    end
else
    if strcmp(aperture, 'coded')
        filePrefix = 'Coded_Orig';
    elseif strcmp(aperture, 'conventional')
        filePrefix = 'Conventional_Orig';
    end
end

%% import, rescale, and crop images 
% output: .mat files
codAp_import_crop_rescale(use_undistorted_images, aperture, numBlurred, filePrefix);


%% generate blur kernels
% codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% larger the window and blur kernel the better!
nums = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
bdims = [33, 27, 33, 33, 37, 51, 49, 51, 57, 61];
KERS = cell(1, 10);
times = zeros(1, 10); 
dim = 500; % dimension of sharp and blurred images

for i = 7
    num = nums(i);
    bdim = bdims(i);
    [x, y] = size(resizedImageData{num});
    midx = round(x/2); midy = round(y/2);

    B = resizedImageData{num}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    S = resizedImageData{1}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
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

%% crop blur kernels
mKers = cell(10, 1);
tempobj = load('mKer_2_1.mat');
mKers{1} = tempobj.mKer_2_1;
tempobj = load('mKer_2_2.mat');
mKers{2} = tempobj.mKer_2_2;
tempobj = load('mKer_2_3.mat');
mKers{3} = tempobj.mKer_2_3;
tempobj = load('mKer_2_4.mat');
mKers{4} = tempobj.mKer_2_4;
tempobj = load('mKer_2_5.mat');
mKers{5} = tempobj.mKer_2_5;
tempobj = load('mKer_2_6.mat');
mKers{6} = tempobj.mKer_2_6;
tempobj = load('mKer_2_7.mat');
mKers{7} = tempobj.mKer_2_7;
tempobj = load('mKer_2_8.mat');
mKers{8} = tempobj.mKer_2_8;
tempobj = load('mKer_2_9.mat');
mKers{9} = tempobj.mKer_2_9;
tempobj = load('mKer_2_10.mat');
mKers{10} = tempobj.mKer_2_10;

pixels = [9, 4, 5, 2, 5, 1, 4, 1, 2, 5]; % approve

for i = 1:length(mKers)
    mKers{i} = cropKernel(mKers{i}, pixels(i));
    figure;
    imshow(mKers{i}./max(max(mKers{i})));
end
save('Coded_Undist_Blur_Kernels.mat', 'mKers');

%% deconvolve, local energy estimate, and depth Map #1 (local independent window)
load('Coded_Undist_Blur_Kernels.mat');
I = imread(strcat(strcat(strcat(lower(filePrefix(1:3)), 'Ap_scene'), num2str(scene)), '_undist.jpg')); 
I_cropped = I(600:1700, 500:2400, :); 
% I_cropped = I(1200:1400, 1200:1400, :);
% I_cropped = I(800:1000, 700:1000, :);
I_gray_double = im2double(rgb2gray(I_cropped));
[ry, cy] = size(I_gray_double);

% deconvSps parameters
we = 0.001;
max_it = 200;

allFocusImage = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
depths = 2.1:0.1:3;
depthMap = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
localEnergies = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))), numBlurred);
for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
    for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
        % re-initialize
        minLocalEnergy = inf;
        minInd = 0;
        minDeconv = zeros(depthMapWindow);
        
        % obtain local window of scene
        I_local = I_gray_double(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1));
        for i = 1:numBlurred
            % rotate kernel to match convention of deconSps %%%
            % tempDeconv = deconvSps(I_local, rot90(mKers{i}, 2), we, max_it);
            tempDeconv = deconvL2(I_local, mKers{i}, we, max_it); % experient with rotation % Use the right orientation Kernel!!!!!!
            % tempDeconv = deconvL2(I_local, rot90(mKers{i}, 2), we);
            
            % calculate reconstruction error for local window
            reconError = I_local - conv2(tempDeconv, mKers{i}, 'same');
            
            % calculate average local energy for this window (Frobenius
            % norm)
            avgLocalEnergy = sum(sum(reconError.^2));
            localEnergies(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1), i) = avgLocalEnergy;
            
            % find minimum energy
            if avgLocalEnergy < minLocalEnergy
                minLocalEnergy = avgLocalEnergy; % update
                minInd = i;
                minDeconv = tempDeconv;
            end
        end
        % update depthMap 
        depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(minInd);
        allFocusImage(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = minDeconv;
        
        fprintf('Row %d of %d and Column %d of %d\n', r, (ry - mod(ry, depthMapWindow(1))), c, (cy - mod(cy, depthMapWindow(2))));
        fprintf('Depth: %d\n\n', depths(minInd));
    end
end
% save('codAp_undist_scene1_depthMap_allFI_localWindow1.mat', 'depthMap', 'allFocusImage', 'localEnergies');


%% deconvolve test scene with all blur kernels #2
load('Coded_Undist_Blur_Kernels.mat');
% I = imread(strcat(strcat(strcat(lower(filePrefix(1:3)), 'Ap_scene'), num2str(scene)), '_undist.jpg'));
I = imread('codAp_scene3_gray_undist.bmp');
I_cropped = I(830:1850, 760:2730, 2);

we = 0.001;
max_it = 200;

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
    deconvImages{i} = deconvL2(im2double(I_cropped), mKers{i}, we, max_it); % Use the right orientation Kernel!!!!!!
    times(i) = toc;
    fprintf('Image %d\n', i);
end


%% local Energy Estimate and depth map #2
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
load('Coded_Undist_Blur_Kernels.mat');
load('codAp_undist_deconv_scene1_deconvL2.mat');
% load('codAp_v2_undist_deconv_scene3_deconvL2.mat');
y_orig = imread(strcat(strcat(strcat(lower(filePrefix(1:3)), 'Ap_scene'), num2str(scene)), '_undist.jpg'));
y_orig = y_orig(600:1700, 500:2400, :);
y = im2double(rgb2gray(y_orig));
% y_orig = imread('codAp_scene3_gray_undist.bmp');
% y_orig = y_orig(830:1850, 760:2730, 2);
% y = im2double(y_orig);
[ry, cy] = size(y);

range = 1:numBlurred; % range of blur kernels using
weights = linspace(1, 1, numBlurred);
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, mKers{i}, 'same');
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
depths = 2.1:0.1:3;
depthMap = zeros((ry - mod(ry, depthMapWindow(1))), (cy - mod(cy, depthMapWindow(2))));
for r = 1:depthMapWindow(1):(ry - mod(ry, depthMapWindow(1)))
    for c = 1:depthMapWindow(2):(cy - mod(cy, depthMapWindow(2)))
        minEng = inf;
        le_plot = zeros(1, length(range));
        for i = range
            curEng = localEnergyEst{i}(r, c);
            le_plot(i) = localEnergyEst{i}(r, c);
            if curEng < minEng
                minEng = curEng;
                depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
        end
        plot(depths(range), le_plot);
    end
end

figure; 
imagesc(depthMap);
axis equal;
title('Depth Map');
colorbar;
figure;
imshow(y_orig);

% 1. try #1 with correct kernel orientation
% 2. use sliding local window and resolve depth for each pixel
% 3. Use more textured objects
% 4. de-rotate (bi-cubic interpolant) (rotate and resize one operation)
% 5  use hamming function
% 6. check linear size of h
% 7. make h much largerrrr
% 8. two lasers for fronto-parallel alignment
% decision tree/forrest 
% compare jpeg with .raw images (.CR2)
% fmincon vs machine learning for local energy weightings


%% generate all-focused image



%% plot results
load('Coded_Undist_Blur_Kernels.mat');
for i = 1:10
   figure;
   imshow(mKers{i}./max(max(mKers{i})));
end

