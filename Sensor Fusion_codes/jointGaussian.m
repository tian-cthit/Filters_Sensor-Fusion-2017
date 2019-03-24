function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
%jointGaussian calculates the joint Gaussian density as defined
%in problem 1.3a. 
%
%Input
%   MU_X        Expected value of x
%   SIGMA2_X    Covariance of x
%   SIGMA2_R    Covariance of the noise r
%
%Output
%   MU          Mean of joint density 
%   SIGMA       Covariance of joint density


%Your code here
A = [1 0;1 1];
b = 0;
mu_xr = [mu_x;0];
sigma_xr = [sigma2_x 0;0 sigma2_r];
% [mu, Sigma] = affineGaussianTransform(mu_xr, sigma_xr, A, b);
mu = A*mu_xr+b;
Sigma = A*sigma_xr*A';
end