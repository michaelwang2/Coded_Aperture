function mKer_cropped = cropKernel(mKer, pix)
% crop border of kernel by pix pixels
mKer_cropped = mKer((1+pix):(end-pix), (1+pix):(end-pix));
end