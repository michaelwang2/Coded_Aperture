function codAp_generate_blur_kernels(blurKerDims, numBlurred, filePrefix)
% this function loads the right cropped images and generate the appropriate
% blur kernels

load(strcat(filePrefix, '_imported_cropped_rescaled.mat'));

blurKernels = cell(numBlurred, 1);
for i = 2:(numBlurred+1)
    blurKernels{i-1} = calcKer(resizedImageData{i}, resizedImageData{1}, blurKerDims);
end

save(strcat(filePrefix, '_blur_kernels'), 'blurKernels');
end