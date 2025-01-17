clc;
clear all

load balloons_ms.mat

p=0.95;

X=normalized(double(Omsi));

Xsize = size(X);
Omega    = (rand(Xsize) >=p); % SR = 0.1
[n1, n2 ,n3]=size(X);
Y = X;

Y2 = Frontal2Lateral(Y); % each lateral slice is a channel of the image
Omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
Omega2(Omega) = Iones(Omega);
Omega2 = Frontal2Lateral(Omega2);
Omega2 = find(Omega2==1);

R = 31;  % sould tune the rank, for CAVE HSI the tuned parameters are : p=0.9 R=31,p=0.95 R=31,p=0.97, R=6,
r = [256 256 3];
opt.maxIter = 1000;
opt.max_mu = 1e10;

opt.rho = 1.01; % should tune the ADMM parameters, here is tuned parameters for CAVE HSI.

opt.mu = 1e-2;


opt.tol = 1e-6;

smooth_para = [0.01 0.01 0.01]; % should tune the smooth regularization paraemeters, here is tuned parameters for CAVE HSI.

[Xhat A B] = old_FSmoothTSVD1(Y2,Omega2,R,[0.00 0.00 0],smooth_para,1,r,opt);
Xhat = Lateral2Frontal(Xhat);
[psnr2, ssim2, fsim2, ergas2, msam2] =MSIQA(Xhat*255,X*255);


