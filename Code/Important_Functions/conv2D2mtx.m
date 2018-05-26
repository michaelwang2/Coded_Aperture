function D = conv2D2mtx(szKer, mask)
% generate matrix form of convolution of mask * Ker
% Ker dimensions must be odd
[r, c] = size(mask);
sizeC = max([szKer(1)-max(0,r-1),szKer(2)-max(0,c-1)],0);
col = zeros(sizeC(1)*sizeC(2), 1); col(1) = mask(1, 1);

start = 1;
D = zeros(sizeC(1)*sizeC(2), szKer(1)*szKer(2));
for i = 1:sizeC(1):(sizeC(1)*sizeC(2))
    temp = zeros(1, szKer(1)*szKer(2));
    count = 1;
   for j = start:szKer(1):(start + szKer(1)*(c-1))
       temp(j:(j+r-1)) = mask(:, count)';
       count = count + 1;
   end
   start = start + szKer(1);
   D(i:(i+sizeC(1)-1), :) = toeplitz(col(i:(i+sizeC(1)-1)), temp);
end
end