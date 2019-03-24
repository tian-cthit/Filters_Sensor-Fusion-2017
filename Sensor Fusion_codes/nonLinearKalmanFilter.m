function [xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, T, Q, S, h, R, type)
%NONLINEARKALMANFILTER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%   T                   Sampling time
%   Q           [n x n] Process noise covariance
%   S           [n x N] Sensor position vector sequence
%   h                   Measurement model function handle
%   R           [n x n] Measurement noise covariance
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%

% Your code here. If you have good code for the Kalman filter, you should re-use it here as
% much as possible.
n = length(x_0);
N = size(Y,2);

xf = zeros(n, N);
Pf = zeros(n,n,N);
xp = zeros(n, N);
Pp = zeros(n,n,N);

[xp(:,1), Pp(:,:,1)] = nonLinKFprediction(x_0, P_0, f, T, Q, type);
[xf(:,1), Pf(:,:,1)] = nonLinKFupdate(xp(:,1), Pp(:,:,1), Y(:,1), S(:,1), h, R, type);

    for i = 2:1:N
        [xp(:,i), Pp(:,:,i)] = nonLinKFprediction(xf(:,i-1), Pf(:,:,i-1), f, T, Q, type);
        [xf(:,i), Pf(:,:,i)] = nonLinKFupdate(xp(:,i), Pp(:,:,i), Y(:,i), S(:,i), h, R, type);
    end




end