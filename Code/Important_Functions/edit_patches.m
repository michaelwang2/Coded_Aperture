%% re-orient, crop, and resize patches

% collect folder names
% rotate and crop images
% folder = dir('*.CR2');
% Patches = cell(length(folder), 1);
load('codAp_v5_R_bk.mat');
figure;
for i = 1:10
    % I = im2double(imread(folder(i).name));
    I = KERS{i};
    fprintf('Crop Image\n');
    % Patches{i} = rgb2gray(imcrop(I_rot));
    max_el = max(max(I));
    I_cropped = max_el.*imcrop(I./max_el);
    Patches{i} = I_cropped./sum(sum(I_cropped));
end

%% resize patches
load('codAp_v5_NR_patterns.mat');
minS = size(Patches{end});

for j = 1:length(Patches)
    Patches{j} = imresize(Patches{j}, minS);
end

save('codAp_v5_NR_patterns_resize.mat', 'Patches');

%% display patches
clear; close; clc;

load('codAp_v5_bk_R_cropped_centered.mat');
for i = 1:length(KERS)
   figure;
   imshow(KERS{i}./max(max(KERS{i})));
end