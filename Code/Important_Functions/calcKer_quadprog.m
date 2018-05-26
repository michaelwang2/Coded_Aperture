%inputs: B,S - gray level blurred and sharp images respectively (double)
%        szKer - 2 element vector specifying the size of the required kernel
%outputs: mKer - the recovered kernel, 
%
%example usage:  mKer = calcKer_quadprog(B, G, [11 11]);

function mKer = calcKer_quadprog(B, S, szKer, szFilt, l1, l2, l3, l4, l5)

% parabolic regularization matrix
R = flipud(vec(parabola_reg(szKer, 2)))';

% Gaussian derivative (sobel) kernels
Gx = conv2(fspecial('gaussian', szFilt, .15), rot90(fspecial('sobel')', 2), 'valid');
Gy = conv2(fspecial('gaussian', szFilt, .15), rot90(fspecial('sobel'), 2), 'valid');
I = zeros(szKer); I(ceil(szKer(1)/2), ceil(szKer(1)/2)) = 1; % Impulse kernel

% re-adjust Blurred image
B = conv2(B, I, 'valid');

% apply Gaussian derivative kernels on images
B_Gx = conv2(B, rot90(Gx, 2), 'valid');
B_Gy = conv2(B, rot90(Gy, 2), 'valid');
S_Gx = conv2(S, rot90(Gx, 2), 'valid');
S_Gy = conv2(S, rot90(Gy, 2), 'valid');

% generate derivative matrix
D = [fliplr(conv2D2mtx(szKer, fspecial('sobel')')); fliplr(conv2D2mtx(szKer, fspecial('sobel')))];

% set up quadratic program for minimization
% 0.5*x'*H*x + f'*x
H = l1.*(im2col(S, szKer, 'sliding')*(im2col(S, szKer, 'sliding')'))...
    + l2.*(im2col(S_Gx, szKer, 'sliding')*(im2col(S_Gx, szKer, 'sliding')'))...
    + l3.*(im2col(S_Gy, szKer, 'sliding')*(im2col(S_Gy, szKer, 'sliding')'))...
    + l4.*((D')*D)...
    + l5.*((R')*R);
f = -(l1.*(im2col(S, szKer, 'sliding')*B(:))...
    + l2.*(im2col(S_Gx, szKer, 'sliding')*B_Gx(:))...
    + l3.*(im2col(S_Gy, szKer, 'sliding')*B_Gy(:)));
A = []; b = [];
Aeq = ones(1, szKer(1)*szKer(2)); beq = 1;
lb = zeros(szKer(1)*szKer(2), 1); ub = ones(szKer(1)*szKer(2), 1);

% solve quadratic program
X = quadprog(H,f,A,b,Aeq,beq,lb,ub);
mKer = reshape(flipud(X), szKer(1), szKer(2));
end

function R = parabola_reg(szKer, pow)
% generate exponential regularization matrix with dimensions dims
R = zeros(szKer(1), szKer(2));
midx = ceil(szKer(1)/2);
midy = ceil(szKer(2)/2);

for i = 1:szKer(1)
    for j = 1:szKer(2)
        R(i, j) = sqrt((i - midy)^2 + (j - midx)^2)^pow;
    end
end
end
