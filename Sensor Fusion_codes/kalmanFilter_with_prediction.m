function [X, P, Xp, Pp] = kalmanFilter_with_prediction(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter. 
%
%Input:
%   Y           [m x N] Measurement sequence
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   A           [n x n] State transition matrix
%   Q           [n x n] Process noise covariance
%   H           [m x n] Measurement model matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   x           [n x N] Estimated state vector sequence
%   P           [n x n x N] Filter error convariance
%

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);

xp = zeros(n,N);
Pp = zeros(n,n,N);
%% 
% calculate the first state
[xp(:,1), Pp(:,:,1)] = linearPrediction(x_0, P_0, A, Q);
[x(:,1), P(:,:,1)] = linearUpdate(xp(:,1), Pp(:,:,1), Y(:,1), H, R);

% state 2 to N
for i = 2:N
    [xp(:,i), Pp(:,:,i)] = linearPrediction(x(:,i-1), P(:,:,i-1), A, Q);
    [x(:,i), P(:,:,i)] = linearUpdate(xp(:,i), Pp(:,:,i), Y(:,i), H, R);
end
X = x;
Xp = xp;
end

function [x, P] = linearPrediction(x, P, A, Q)
x = A*x;
P = A*P*A' + Q;
end


function [x, P] = linearUpdate(x, P, y, H, R)
S = H*P*H' + R;
v = y - H*x;
K = P*H'/S;

x = x + K*v;
P = P - K*S*K';
end