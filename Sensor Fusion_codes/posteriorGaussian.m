function [mu, sigma2] = posteriorGaussian(mu_p, sigma2_p, y, sigma2_r)
%posteriorGaussian performs a single scalar measurement update with a
%measurement model which is simply "y = x + noise".
%
%Input
%   MU_P            The mean of the (Gaussian) prior density.
%   SIGMA2_P        The variance of the (Gaussian) prior density.
%   SIGMA2_R        The variance of the measurement noise.
%   Y               The given measurement.
%
%Output
%   MU              The mean of the (Gaussian) posterior distribution
%   SIGMA2          The variance of the (Gaussian) posterior distribution

%Your code here
A = [1 0;1 1];
b = 0;
mu_xr = [mu_p;0];
sigma_xr = [sigma2_p 0;0 sigma2_r];
mu_xy = A*mu_xr+b;
Sigma_xy = A*sigma_xr*A';

mu = mu_p + Sigma_xy(1,1)*Sigma_xy(2,2)^-1*(y - mu_xy(2));
sigma2 = Sigma_xy(1,1)-Sigma_xy(1,2)*Sigma_xy(2,2)^-1*Sigma_xy(2,1);

end