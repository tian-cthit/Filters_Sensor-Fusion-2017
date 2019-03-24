clear
clc

x = 20; %assume mean of real temperature to be 20
mu = [x;0];

Sigma = [19.5^2 0;0 2^2];
%% a
[mu_XY, Sigma_XY] = jointGaussian(x, 19.5^2, 2^2); % joint density
[ xy ,p ] = sigmaEllipse2D( mu_XY, Sigma_XY);       %3-sigma ellipse
%% b
y1 = 17;
y2 = 27;

[mu_1, sigma_1] = posteriorGaussian(x, 19.5^2, y1, 2^2);        %posterior distribution
[mu_2, sigma_2] = posteriorGaussian(x, 19.5^2, y2, 2^2);
plot(-10:0.05:50,normpdf(-10:0.05:50,mu_1,sigma_1),'r');        %plot  posterior distribution
hold on;
plot(-10:0.05:50,normpdf(-10:0.05:50,mu_2,sigma_1),'b');
legend('y1=17','y2=27');