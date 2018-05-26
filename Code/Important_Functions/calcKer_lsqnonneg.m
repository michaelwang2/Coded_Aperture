%inputs: B,S - gray level blurred and sharp images respectively (double)
%        szKer - 2 element vector specifying the size of the required kernel
%outputs: mKer - the recovered kernel, 

%example usage:  mKer = calcKer_lsqnonneg(B, S, [11 11]);

function mKer = calcKer_lsqnonneg(B, S, szKer, l1, l2)
%get the "valid" pixels from B (i.e. those that do not depend 
%on zero-padding or a circular assumption
imBvalid = B(ceil(szKer(1)/2):end-floor(szKer(1)/2), ...
  ceil(szKer(2)/2):end-floor(szKer(2)/2));

%get a matrix where each row corresponds to a block from S 
%the size of the kernel
S_Matrix = im2col(l1.*S, szKer, 'sliding')';

% generate derivative matrix
D = fliplr(conv2D2mtx(szKer, l2.*fspecial('log')'));
bzeros = zeros(size(D, 1), 1);

% solve using non-negative least squares
C = [S_Matrix;...
     D;...
     ones(1, size(D, 2))];
d = [l1.*imBvalid(:);...
    bzeros;...
    1];
vXcorrKer = lsqnonneg(C, d);

%reshape and rotate 180 degrees to get the convolution kernel (preflip)
mKer = rot90(reshape(vXcorrKer, szKer), 2);
end