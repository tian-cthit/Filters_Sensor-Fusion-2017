function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f



if nargin < 4
    N = 5000;
end

%Your code here
R = mvnrnd(mu_x,Sigma_x,N); 
Sample = cell(1,N);
for i = 1:1:N
    Sample{i} = f(R(i,:)');
end
ESum = zeros(length(Sample{1}),1);
CSum = zeros(length(Sample{1}),length(Sample{1}));
for i = 1:1:N
    ESum = ESum + Sample{i};
end
E = ESum/N;
for i = 1:1:N
    CSum = CSum + (Sample{i} - E)*(Sample{i} - E)';
end
C = CSum/N;
mu_y = E;
Sigma_y = C;
y = zeros(N,2);
for i = 1:1:N
y(i,:) = Sample{i};
end
y_s = y;

end