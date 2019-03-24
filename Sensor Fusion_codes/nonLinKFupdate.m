function [x, P] = nonLinKFupdate(x, P, y, s, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   s           [2 x 1] sensor position vector
%   h           Measurement model function handle
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

n = length(x);
m = length(y);

switch type
    case 'EKF'
        [y_2, C] = h(x, s);
        % Your EKF update here
        S = C*P*C' + R;
        K = P*C'*S^-1;
        x = x + K*(y-y_2);
        P = P - K*S*K';
    case 'UKF'
        [SP,W] = sigmaPoints(x, P, type);
        
        y_k = zeros(m, 1);
        P_xy = zeros(n, m);
        S = zeros(m,m);
        
        for i=1:1:2*n+1
            [y_2, C] = h(SP(:,i), s);
            y_k = y_k + y_2*W(i);
        end
        
        for i=1:1:2*n+1
            [y_2, C] = h(SP(:,i), s);
            P_xy = P_xy + (SP(:,i) - x)*(y_2 - y_k)'*W(i);
            S = S + (y_2 - y_k)*(y_2 - y_k)'*W(i);
        end
        
        S = S + R;
        x = x + P_xy*(S^-1)*(y - y_k);
        P = P - P_xy*(S^-1)*P_xy';
        
    case 'CKF'
        [SP,W] = sigmaPoints(x, P, type);
        
        y_k = zeros(m, 1);
        P_xy = zeros(n, m);
        S = zeros(m,m);

        for i=1:1:2*n
            [y_2, C] = h(SP(:,i), s);
            y_k = y_k + y_2*W(i);
        end
        
        for i=1:1:2*n
            [y_2, C] = h(SP(:,i), s);
            P_xy = P_xy + (SP(:,i) - x)*(y_2 - y_k)'*W(i);
            S = S + (y_2 - y_k)*(y_2 - y_k)'*W(i);
        end
        
        S = S + R;
        x = x + P_xy*(S^-1)*(y - y_k);
        P = P - P_xy*(S^-1)*P_xy';
        % Your CKF update here
        
    otherwise
        error('Incorrect type of non-linear Kalman filter')
end

end
