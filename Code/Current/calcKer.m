%inputs: B,G - gray level blurred and sharp images respectively (double)
%        szKer - 2 element vector specifying the size of the required kernel
%outputs: mKer - the recovered kernel, 
%         imBsynth - the sharp image convolved with the recovered kernel
%
%example usage:  mKer = calcKer(B, G, [11 11]);

function [mKer, imBsynth] = calcKer(B, G, szKer)

  %get the "valid" pixels from B (i.e. those that do not depend 
  %on zero-padding or a circular assumption
  imBvalid = B(ceil(szKer(1)/2):end-floor(szKer(1)/2), ...
      ceil(szKer(2)/2):end-floor(szKer(2)/2));

  %get a matrix where each row corresponds to a block from G 
  %the size of the kernel
  mGconv = im2col(G, szKer, 'sliding')';

  %solve the over-constrained system using MATLAB's backslash
  %to get a vector version of the cross-correlation kernel
  vXcorrKer = mGconv \ imBvalid(:);

  %reshape and rotate 180 degrees to get the convolution kernel
  mKer = rot90(reshape(vXcorrKer, szKer), 2);

  if (nargout > 1)
      %if there is indeed a convolution relationship between B and G
      %the following will result in an image similar to B
      imBsynth = conv2(G, mKer, 'valid');
  end

end