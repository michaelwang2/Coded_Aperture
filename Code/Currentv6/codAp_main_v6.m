%% generate blur kernels
% codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% larger the window and blur kernel the better!
clear; close; clc;

load('Resizing/codAp_v5_R_patterns_resize.mat');

nums = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
bdims = [27 27 33 33 37 45 47 47 53 53];
KERS = cell(1, 10);
times = zeros(1, 10); 
dim = 500; % dimension of sharp and blurred images

for i = 10
    num = nums(i);
    bdim = bdims(i);
    [x, y] = size(Patches{num});
    midx = round(x/2); midy = round(y/2);

    B = Patches{num}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
    S = Patches{1}((midx-dim/2):(midx+dim/2 - 1), (midy-dim/2):(midy+dim/2 - 1));
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
    l1 = 1;
    l2 = 0.1;
    mKer = calcKer_lsqnonneg(B, S, szKer, l1, l2);
    times(i) = toc;
    KERS{i} = mKer;
    fprintf('Kernel number %d is done!\n', i);
end

%% deconvolve test scene 1
clear; close; clc;
load('Kernels/codAp_v5_bk_NR_cropped.mat');
I = imread('IMG_0140.CR2');
I_cropped = I(600:1550, 800:2400, :);

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

%% local Energy Estimate and depth map sliding window
clear; close; clc;

% scene number
scene = 1;

load('codAp_v6_1_3_bks.mat');
load(['codAp_v6_1_3_deconvL2_flipped_scene', num2str(scene), '.mat']);
load('1_3_focus/codAp_v6_1_3_scenes.mat');

numBlurred = length(KERS);
depthMapWindow = [61,61];  
reconErrors = cell(numBlurred, 1);
localEnergyEst = cell(numBlurred, 1);
SSIMest = cell(numBlurred, 1);
y = scenes{scene};
[ry, cy] = size(y);

range = 4:10;
y_other = cell(numBlurred, 1);
for i = range
    y_other{i} = conv2(deconvImages{i}, KERS{i}, 'same');
    reconErrors{i} = y - y_other{i};
end

% generate depth map
figure; hold on;
% depths = (2.1:0.1:3); % m 
depths = 1.4:0.1:2.3;
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
depthMapSSIM = zeros(ry, cy);
K1 = 1e-3; K2 = 1e-3;
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
imshow(y);
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

%% Simulate PSF Diameter
clear; close; clc;

% toggles
focus_plane_lock = 2000; % mm
od = linspace(focus_plane_lock, 20000, 100); % object distance mm
percentage = 0.85;

% calculate distance between lens and image sensor
% f = 50; % mm
f = 80;
s = 1/((1/f) - (1/focus_plane_lock)); % mm

% calculate PSF diameter
pixel_width = 14.8/2304; % mm/pixel
D = 11; % aperture diameter mm
v = 1./((1/f) - (1./od)); % image distance mm
dm = D.*abs(s - v)./v; % PSF diameter mm
dp = dm./pixel_width; % pixels
max_PSF = (D*abs(s - f)/f)/pixel_width; % max PSF diameter in pixels
dD_da = @(a, f, fp, Dap)(Dap.*((a.*f)./((a + fp).^2 .*(f - fp)) - f./((a + fp).*(f - fp))));

% plot
figure; hold on;
plot((od-focus_plane_lock)./1000, dp, 'k-'); 
plot([0, od(end)-focus_plane_lock]./1000, [max_PSF, max_PSF], 'r--');
plot([0, od(end)-focus_plane_lock]./1000, [max_PSF, max_PSF].*percentage, 'b--');
plot([focus_plane_lock, focus_plane_lock]./1000 + 1, [0, max(dp)], 'g--');
grid on; box on; 
xlabel('Distance from Focal Plane [m]');
ylabel('PSF Diameter [pixels]');
title(['Focus locked at ', num2str(focus_plane_lock/1000), ' m']);
legend('PSF Diameter', 'Theoretical Max', [num2str(percentage*100), '\% of Theoretical Max'], 'Location', 'SouthEast');

% plot
figure;
deriv = dD_da(od-focus_plane_lock, f, focus_plane_lock, D);
plot((od-focus_plane_lock)./1000, deriv, 'k-');
grid on; box on; 
xlabel('Distance from Focal Plane [m]');
yy=ylabel('$$\frac{d D_{PSF}}{d(o_d - f_p)}$$', 'rotation', 0);
title(['Focus locked at ', num2str(focus_plane_lock/1000), ' m']);
set(yy, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

%% Simulate PSF Diameter with focus plane lock sweep
clear; close; clc;

% toggles
focus_plane_locks = linspace(1000, 3000, 11); % mm

figure; hold on;
for i = 1:length(focus_plane_locks)
    focus_plane_lock = focus_plane_locks(i);
    od = linspace(focus_plane_lock, focus_plane_lock+4000, 100); % object distance mm

    % calculate distance between lens and image sensor
    % f = 50; % mm
    f = 80;
    s = 1/((1/f) - (1/focus_plane_lock)); % mm

    % calculate PSF diameter
    pixel_width = 14.8/2304; % mm/pixel
    D = 11; % aperture diameter mm
    v = 1./((1/f) - (1./od)); % image distance mm
    dm = D.*abs(s - v)./v; % PSF diameter mm
    dp = dm./pixel_width; % pixels
    max_PSF = (D*abs(s - f)/f)/pixel_width; % max PSF diameter in pixels

    % plot
    plot((od-focus_plane_lock)./1000, dp);
end
grid on; box on; 
xlabel('Distance from Focal Plane [m]');
ylabel('PSF Diameter [pixels]');
title('Focus Distance Lock Sweep');
legend('1 m', '1.2 m', '1.4 m', '1.6 m', '1.8 m', '2.0 m', '2.2 m', '2.4 m'...
    , '2.6 m', '2.8 m', '3.0 m');

%%
clear; close all; clc;

load('codAp_v7_bks_cropped.mat');
n = size(KERS{end}, 1);
template = zeros(n, n);
paddedKers = zeros(n,n,10);
for i = 1:9
    tt = template;
    m = size(KERS{i}, 1);
    wing = (n-m)/2;
    tt((wing+1):(n-wing), (wing+1):(n-wing)) = KERS{i}./max(max(KERS{i}));;
    paddedKers(:, :, i) = tt;
end

figure;
montage(paddedKers, 'size', [2,5]);
