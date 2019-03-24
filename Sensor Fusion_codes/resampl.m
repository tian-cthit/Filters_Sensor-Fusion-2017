function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Your code here!
N = length(W);
n = size(X,1);

u = ([0:N-1]+rand(1))/N;
W_batman = cumsum(W);       %interval weights
W_batman = W_batman/W_batman(N);   %normalization

[dum,ind1] = sort([u,W_batman]);
ind2 = find(ind1<=N);
j = ind2 - (0:N-1);   %index
Xr = X(:,j);

Wr = ones(1,N)./N;
end
