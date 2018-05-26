function S = SSIM(O, B, K1, K2)
% finds the structure similarity index between Original Image (O) and
% Blurred image (B)

% average luminance
mu_o = sum(sum(O))./numel(O);
mu_b = sum(sum(B))./numel(B);

% standard deviations
C = cov(O(:), B(:)); % covariance matrix
sig_o = sqrt(C(1,1));
sig_b = sqrt(C(2,2));
sig_ob = C(1, 2); 


% constants
L = 1; % dynamics range assuming double
C1 = (K1*L)^2;
C2 = (K2*L)^2;

% SSIM 
S = (2*mu_o*mu_b + C1)*(2*sig_ob + C2)/((mu_o^2 + mu_b^2 + C1)*(sig_o^2 + sig_b^2 + C2));
end