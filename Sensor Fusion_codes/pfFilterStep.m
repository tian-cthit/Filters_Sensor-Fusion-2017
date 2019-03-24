function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, y_k, proc_f, proc_Q, meas_h, meas_R)
%PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
%
% Input:
%   X_kmin1     [n x N] Particles for state x in time k-1
%   W_kmin1     [1 x N] Weights for state x in time k-1
%   y_k         [m x 1] Measurement vector for time k
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%
% Output:
%   X_k         [n x N] Particles for state x in time k
%   W_k         [1 x N] Weights for state x in time k

% Y-our code here!
n = size(X_kmin1,1);
N = size(X_kmin1,2);
X_k = zeros(n,N);
Y_k_pre = zeros(length(y_k),1);
for i  = 1:1:N
    X_k(:,i) = proc_f(X_kmin1(:,i));
    X_k(:,i) = mvnrnd(X_k(:,i),proc_Q);
    Y_k_pre(:,i) = meas_h(X_k(:,i));    %m x N    
end

Y = mvnpdf(y_k',Y_k_pre',meas_R);  % Y:Nx1
W_k_0 = W_kmin1 .* Y';   
W_k_sum = sum(W_k_0);
W_k = W_k_0/W_k_sum;
end