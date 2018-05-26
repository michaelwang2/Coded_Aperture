% convert .CR2 images to gray (.png)

dirname = 'RAW/';
imageNames = dir(strcat(dirname, '*.CR2'));

for i = 1:length(imageNames)
    raw = imread(strcat(dirname, imageNames(i).name));
    gray = rgb2gray(raw);
    imwrite(gray, strcat(dirname, imageNames(i).name(1:(end-4)), '_gray', '.bmp'), 'bmp');
end