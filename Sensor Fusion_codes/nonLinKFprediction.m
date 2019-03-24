function [x, P] = nonLinKFprediction(x, P, f, T, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%   T           Sampling time
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
n = length(x);
switch type
    case 'EKF'
        
        % Your EKF code here
        [x, A] = f(x,T);
        P = A*P*A' + Q;
        
    case 'UKF'
        % Your UKF code here
        [SP,W] = sigmaPoints(x, P, type);
        
        x = zeros(n,1);
        for i = 1:1:2*n+1
            x = x + f(SP(:,i),T)*W(i);
        end
            
        P = zeros(n,n);
        for i = 1:1:2*n+1
            P = P + (f(SP(:,i),T) - x)*(f(SP(:,i),T) - x)'*W(i);
        end
            P = P + Q;
    case 'CKF'
        
        % Your CKF code here
        [SP,W] = sigmaPoints(x, P, type);
        
        x = zeros(n,1);
        for i = 1:1:2*n
            x = x + f(SP(:,i),T)*W(i);
        end
            
        P = zeros(n,n);
        for i = 1:1:2*n
            P = P + (f(SP(:,i),T) - x)*(f(SP(:,i),T) - x)'*W(i);
        end
            P = P + Q;
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end
end