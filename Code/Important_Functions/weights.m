%% Create energy dataset 

clear; close; clc;

addpath('Calibration/Weight');

load('Patches_rearrange.mat');
load('Kernels/codAp_v4_bk_cropped.mat');

range = 2:10; % range of kernels to use
window = [19, 19];
kernels = cell(length(range), 1);
count = 1;
for i = range
    kernels{count} = KERS{i};
    count = count + 1;
end

% depths
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m

% form energies matrix
we = 0.006;
max_it = 200;
spacing = 2; % spacing for sliding window
p = length(Patches_rearrange); % number of patches to use
correctDepths = [];
localEnergies = []; 
for i = 1:p
    fprintf('Patch Number %i\n', i);
    for k = range + 1
        fprintf('Depth Number %i\n', k);
        energies = bk2energy(Patches_rearrange{i}{k}, kernels, window, we, max_it, spacing);
        localEnergies = [localEnergies; energies];
        correctDepths = [correctDepths; ones(size(energies, 1), 1).*depths(k-1)];
    end
    fprintf('\n');
end

save('weightings9.mat', 'localEnergies', 'correctDepths', 'depths', 'range', 'window');

%% Create energy dataset (flipped kernels for deconvL2)

clear; close; clc;

addpath('Calibration/Weight');

load('Patches_rearrange.mat');
load('Kernels/codAp_v4_bk_cropped.mat');

range = 2:10; % range of kernels to use
window = [19, 19];
kernels = cell(length(range), 1);
count = 1;
for i = range
    kernels{count} = KERS{i};
    count = count + 1;
end

% depths
lens_length_correction = -0.0209;% m
depths = (2.1:0.1:3) + lens_length_correction; % m

% form energies matrix
we = 0.006;
max_it = 200;
spacing = window(1); % spacing for sliding window
p = length(Patches_rearrange); % number of patches to use
correctDepths = [];
localEnergies = []; 
for i = 1:p
    fprintf('Patch Number %i\n', i);
    for k = range + 1
        fprintf('Depth Number %i\n', k);
        energies = bk2energy_flipped(Patches_rearrange{i}{k}, kernels, window, we, max_it, spacing);
        localEnergies = [localEnergies; energies];
        correctDepths = [correctDepths; ones(size(energies, 1), 1).*depths(k-1)];
    end
    fprintf('\n');
end

save('weightings8.mat', 'localEnergies', 'correctDepths', 'depths', 'range', 'window');

%% Train weightings with no bias
clear; close; clc;

addpath('Calibration/Weight');

load('weightings8.mat');

% w0 = ones(length(range), 1); % w0 = linspace(1.5, 1, length(range));

% generate Latin Hyperube Samples
numSamples = 100; 
w0s = lhsdesign(numSamples, length(range));
w0s = w0s.*(ones(1, length(range)).*2.5) + repmat(ones(1, length(range)).*0.5, numSamples, 1);
minVal = inf; 
for i = 1:numSamples
    fprintf('Number %i\n', i);
    w0 = w0s(i, :);
    w = zeros(range(end), 1);
    [ws, zval] = calcWeightings(localEnergies, depths(range), correctDepths, w0);
    w(range) = ws;
    if zval < minVal
        minVal = zval;
        minw = w;
    end
end

save('w8_19x19_deconvL2_2-10_nobias_flipped.mat', 'minw', 'window', 'range');

%% Train weightings with bias
clear; close; clc;

addpath('Calibration/Weight');

load('weightings8.mat');

% w0 = ones(length(range), 1); % w0 = linspace(1.5, 1, length(range))';
% b0 = ones(length(range), 1).*0;

% generate Latin Hyperube Samples
numSamples = 100; 
w0s = lhsdesign(numSamples, length(range));
w0s = w0s.*(ones(1, length(range)).*2.5) + repmat(ones(1, length(range)).*0.5, numSamples, 1);
b0s = lhsdesign(numSamples, length(range));
b0s = b0s.*(ones(1, length(range)).*0.5);
minVal = inf;
for i = 1:numSamples
    fprintf('Number %i\n', i);
    z0 = [w0s(i, :)'; b0s(i, :)'];
    w = zeros(range(end), 1); b = zeros(range(end), 1);
    [ws, bs, zval] = calcWeightings_bias(localEnergies, depths, correctDepths, z0);
    w(range) = ws;
    b(range) = bs;
    if zval < minVal
        minVal = zval;
        minw = w;
        minb = b;
    end
end

save('w8_19x19_deconvL2_2-10_bias_flipped.mat', 'minw', 'minb', 'window', 'range');

% 1. get energy dataset using local window method 
% optimization to choose the best one
% 3. get higher resolution calibration photos
% 4. use deconvSps
% 5. A PSF, which yields a deblurred image with the best quality is
% chosen as the right PSF
% 7. Image Patches: blur them for each kernel and deblur with each kernel
% to get reconstruction error (Dsiadvantage: does not calibrate for test
% scene conditions)
% 8. Use the size of blur kernel as "sliding window" for energy calc
% 9. do better camera focusing
% 10. larger sliding window (25 x 25)
% 11. get blur kernels for sub windows of calibration pattern
% 12. Support Vector Machines
% 13. Use conference room TV as calibration screen for everything
% 14. What should be the fixed focus distance?????? How does focus distance
% affect optimal range of camera?
% 15. Do blur kernels have to take into account object scaling??
% 16. Compute weightings such that it normalizes energy term. 
% 17. Set higher exposure time
% 18. Try PSF estimation using point light source
% 19. When using deconvL2, deconvolve full image and THEN crop. 
% 20. DECONVL2 REQUIRES FLIPPED KERNELS
% 21. Have shutter actually close fully during shot
% 22. Have the central pixel be corners of 4 windows. Choose window with
% best error. 
% 23. segmentation + depthmap
% 24. use SIFT, RIFT to compare windows
% 25. Combining features to make more robust classification
% 26. Histogram of Gradients (HOG)
% 27. Read Image Patch Similarity Paper
% 28. Try local sliding window again but use 'valid' instead


%% Support Vector Machines 3D Visualization
clear; close; clc;

n = 50000;
random_points = rand(n, 3);
[val, inds] = min(random_points, [], 2);
depthsx = random_points(inds == 1, :);
depthsy = random_points(inds == 2, :);
depthsz = random_points(inds == 3, :);

% plot
figure; hold on;
xx = plot3(depthsx(:, 1), depthsx(:, 2), depthsx(:, 3), 'r.');
% yy = plot3(depthsy(:, 1), depthsy(:, 2), depthsy(:, 3), 'g.');
% zz = plot3(depthsz(:, 1), depthsz(:, 2), depthsz(:, 3), 'b.');
grid on; box on; axis equal;
view(45, 45);
title('Classification');
xlabel('Energies X');
ylabel('Energies Y');
zlabel('Energies Z');
legend([xx yy zz], 'Depth X', 'Depth Y', 'Depth Z');

%% Support Vector Machines 3D with 2 considered Visualization
clear; close; clc;

n = 50000;
random_points = rand(n, 3);
[val, inds] = min(random_points(:, 1:2), [], 2);
depthsx = random_points(inds == 1, :);
depthsy = random_points(inds == 2, :);

% plot
figure; hold on;
xx = plot3(depthsx(:, 1), depthsx(:, 2), depthsx(:, 3), 'r.');
yy = plot3(depthsy(:, 1), depthsy(:, 2), depthsy(:, 3), 'g.');
grid on; box on; axis equal;
view(45, 45);
title('Classification');
xlabel('Energies X');
ylabel('Energies Y');
zlabel('Energies Z');
legend([xx yy zz], 'Depth X', 'Depth Y', 'Depth Z');



