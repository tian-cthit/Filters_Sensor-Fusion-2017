function [ xy ,p ] = sigmaEllipse2D( mu, Sigma, level, npoints )
%SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
% a sigma level in the Gaussian density defined by mean and covariance.
%
%Input:
%   MU          [2 x 1] Mean of the Gaussian density
%   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
%   LEVEL       Which sigma level curve to plot. Default = 3.
%   NPOINTS     Number of points on the ellipse to generate. Default = 32.
%
%Output:
%   XY          2 x npoints matrix. First row holds x-coordinates, second
%               row holds the y-coordinates.

if nargin < 3
    level = 3;
end
if nargin < 4
    npoints = 32;
end

%Your code here
i=1;
xy=zeros(2,32);
for theta = 0:2*pi/(npoints-1):2*pi
   xy(:,i) = mu + 3*sqrtm(Sigma)*[cos(theta);sin(theta)];
   i=i+1;
end
p = plot(xy(1,:),xy(2,:));


end
