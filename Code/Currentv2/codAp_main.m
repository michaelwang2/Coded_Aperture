%% import, rescale, and crop images 
% output: .mat files
% [XMIN YMIN WIDTH HEIGHT]
clear; close; clc;
% import
I0 = imread('codAp_2_0_gray_undist.bmp');
I1 = imread('codAp_2_1_gray_undist.bmp');
I2 = imread('codAp_2_2_gray_undist.bmp');
I3 = imread('codAp_2_3_gray_undist.bmp');
I4 = imread('codAp_2_4_gray_undist.bmp');
I5 = imread('codAp_2_5_gray_undist.bmp');
I6 = imread('codAp_2_6_gray_undist.bmp');
I7 = imread('codAp_2_7_gray_undist.bmp');
I8 = imread('codAp_2_8_gray_undist.bmp');
I9 = imread('codAp_2_9_gray_undist.bmp');
I10 = imread('codAp_3_0_gray_undist.bmp');

% crop
codAp_2_0_gray_undist = imcrop(im2double(I0(:, :, 2)), [323, 41, 2847, 2127]);
codAp_2_1_gray_undist = imcrop(im2double(I1(:, :, 2)), [355, 88, 2711, 2023]);
codAp_2_2_gray_undist = imcrop(im2double(I2(:, :, 2)), [537, 128, 2597, 1929]);
codAp_2_3_gray_undist = imcrop(im2double(I3(:, :, 2)), [497, 175, 2486, 1847]);
codAp_2_4_gray_undist = imcrop(im2double(I4(:, :, 2)), [553, 214, 2384, 1767]);
codAp_2_5_gray_undist = imcrop(im2double(I5(:, :, 2)), [630, 225, 2287, 1697]);
codAp_2_6_gray_undist = imcrop(im2double(I6(:, :, 2)), [693, 248, 2202, 1633]);
codAp_2_7_gray_undist = imcrop(im2double(I7(:, :, 2)), [637, 255, 2131, 1581]);
codAp_2_8_gray_undist = imcrop(im2double(I8(:, :, 2)), [589, 268, 2069, 1531]);
codAp_2_9_gray_undist = imcrop(im2double(I9(:, :, 2)), [708, 278, 1993, 1486]);
codAp_3_0_gray_undist = imcrop(im2double(I10(:, :, 2)), [701, 297, 1930, 1433]);

% B = imresize(A, [NUMROWS NUMCOLS])
[r, c] = size(codAp_3_0_gray_undist);
resize_2_0 = imresize(codAp_2_0_gray_undist, [r, c], 'bicubic');
resize_2_1 = imresize(codAp_2_1_gray_undist, [r, c], 'bicubic');
resize_2_2 = imresize(codAp_2_2_gray_undist, [r, c], 'bicubic');
resize_2_3 = imresize(codAp_2_3_gray_undist, [r, c], 'bicubic');
resize_2_4 = imresize(codAp_2_4_gray_undist, [r, c], 'bicubic');
resize_2_5 = imresize(codAp_2_5_gray_undist, [r, c], 'bicubic');
resize_2_6 = imresize(codAp_2_6_gray_undist, [r, c], 'bicubic');
resize_2_7 = imresize(codAp_2_7_gray_undist, [r, c], 'bicubic');
resize_2_8 = imresize(codAp_2_8_gray_undist, [r, c], 'bicubic');
resize_2_9 = imresize(codAp_2_9_gray_undist, [r, c], 'bicubic');
resize_3_0 = codAp_3_0_gray_undist;
Images = {resize_2_0, resize_2_1, resize_2_2, resize_2_3, resize_2_4,...
    resize_2_5, resize_2_6, resize_2_7, resize_2_8, resize_2_9, resize_3_0};

% save
save('codAp_v2_raw_undist.mat', 'Images');


%% generate blur kernels
% codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% larger the window and blur kernel the better!
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

%% deconvolve, local energy estimate, and depth Map #1 (local independent window) v2_v1
clear; close; clc;
load('Kernels/v2_v1_blur_kernels_cropped.mat');
I = imread('RAW/codAp_scene2.CR2');
I_cropped = I(500:1700, 600:2500, :);
I_gray_double = im2double(rgb2gray(I_cropped));
[ry, cy] = size(I_gray_double);

% deconvSps parameters
we = 0.001;
max_it = 200;

dmpsize = 15;
numBlurred = length(v2_v1);
range = 3:numBlurred;
depthMapWindow = [dmpsize, dmpsize];
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
        for i = range
            % rotate kernel to match convention of deconSps %%%
            % tempDeconv = deconvSps(I_local, rot90(mKers{i}, 2), we, max_it);
            tempDeconv = deconvL2(I_local, v2_v1{i}, we, max_it); % experient with rotation % Use the right orientation Kernel!!!!!!
            % tempDeconv = deconvL2(I_local, rot90(mKers{i}, 2), we);
            
            % calculate reconstruction error for local window
            reconError = I_local - conv2(tempDeconv, v2_v1{i}, 'same');
            
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

%% deconvolve, local energy estimate, and depth Map #1 (local independent window) v2_v3
clear; close; clc;
load('Kernels/v2_v3_blur_kernels_cropped.mat');
I = imread('RAW/codAp_scene2.CR2');
I_cropped = I(500:1700, 600:2500, :);
I_gray_double = im2double(rgb2gray(I_cropped));
[ry, cy] = size(I_gray_double);

% deconvSps parameters
we = 0.001;
max_it = 200;

dmpsize = 15;
numBlurred = length(v2_v3);
depthMapWindow = [dmpsize, dmpsize];
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
            tempDeconv = deconvL2(I_local, v2_v3{i}, we, max_it); % experient with rotation % Use the right orientation Kernel!!!!!!
            % tempDeconv = deconvL2(I_local, rot90(mKers{i}, 2), we);
            
            % calculate reconstruction error for local window
            reconError = I_local - conv2(tempDeconv, v2_v3{i}, 'same');
            
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


%% deconvolve test scene with all blur kernels #2 v2_v1
clear; close; clc;
load('Kernels/v2_v1_blur_kernels_cropped.mat');
I = imread('RAW/codAp_scene6.CR2');
I_cropped = I(500:1700, 600:2500, :);

we = 0.001;
max_it = 200;

numBlurred = length(v2_v1);
deconvImages = cell(numBlurred, 1);
times = zeros(numBlurred, 1);
for i = 1:numBlurred
    tic;
    % rotate kernel to match convention of deconSps %%%
    % deconvolve window by window so the artifacts do not affect other
    % scenes
    % conv2 rotates the kernel / filter mask by 180 degrees before performing filtering
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
    deconvImages{i} = deconvL2(im2double(rgb2gray(I_cropped)), v2_v1{i}, we, max_it); % Use the right orientation Kernel!!!!!!
    times(i) = toc;
    fprintf('Image %d\n', i);
end


%% local Energy Estimate and depth map #2 v2_v1
clear; close; clc;
load('Kernels/v2_v1_blur_kernels_cropped.mat');
load('v2_v1_codAp_scene2_deconvImages.mat');

numBlurred = length(v2_v1);
dmpsize = 21;
depthMapWindow = [dmpsize, dmpsize];
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
y_orig = imread('RAW/codAp_scene2.CR2');
y_orig = y_orig(500:1700, 600:2500, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

range = 1:numBlurred; % range of blur kernels using
weights = linspace(1, 1, numBlurred);
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, v2_v1{i}, 'same');
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

% 1. try #1 with correct kernel orientation
% 2. use sliding local window and resolve depth for each pixel
% 6. check linear size of h
% 8. two lasers for fronto-parallel alignment
% decision tree/forrest 
% fmincon vs machine learning for local energy weightings


%% deconvolve test scene with all blur kernels #2 v2_v3
clear; close; clc;
load('Kernels/v2_v3_blur_kernels_cropped.mat');
I = imread('RAW/codAp_scene6.CR2');
I_cropped = I(500:1700, 600:2500, :);

we = 0.001;
max_it = 200;

numBlurred = length(v2_v3);
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
    deconvImages{i} = deconvL2(im2double(rgb2gray(I_cropped)), v2_v3{i}, we, max_it); % Use the right orientation Kernel!!!!!!
    times(i) = toc;
    fprintf('Image %d\n', i);
end


%% local Energy Estimate and depth map #2 v2_v3
clear; close; clc;
load('Kernels/v2_v3_blur_kernels_cropped.mat');
load('v2_v3_codAp_scene2_deconvImages.mat');
load('w4_v3.mat');

numBlurred = length(v2_v3);
dmpsize = 13;
depthMapWindow = [dmpsize, dmpsize];
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
y_orig = imread('RAW/codAp_scene2.CR2');
y_orig = y_orig(500:1700, 600:2500, :);
y = im2double(rgb2gray(y_orig));
[ry, cy] = size(y);

range = 1:numBlurred; % range of blur kernels using
weights = linspace(1, 1, numBlurred);
for i = range
    reconErrors{i} = y - conv2(deconvImages{i}, v2_v3{i}, 'same');
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
            le_plot(i) = localEnergyEst{i}(r, c);
            curEng = localEnergyEst{i}(r, c);
            if curEng < minEng
                minEng = curEng;
                depthMap(r:(r + depthMapWindow(1) - 1), c:(c + depthMapWindow(2) - 1)) = depths(i);
            end
        end
        plot(depths(range), le_plot);
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
