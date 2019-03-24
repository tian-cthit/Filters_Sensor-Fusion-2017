function [xs, Ps, xf, Pf, xp, Pp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, sigmaPoints, type)
%NONLINRTSSMOOTHER Filters measurement sequence Y using a 
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
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%   xs          [n x N]     Smoothed estimates for times 1,...,N
%   Ps          [n x n x N] Smoothing error convariance

% your code here!
% We suggest that you copy the code from your non-linear Kalman filter 
% and then add the smoothing code afterwards. You can also call the 
% non-linear Kalman filter function from here, but then you need to
% modify its input parameter list to include the sigmaPoints function
% handle, and include its function implementation below.'
n = length(x_0);
N = size(Y,2);

xf = zeros(n, N);
Pf = zeros(n,n,N);
xp = zeros(n, N);
Pp = zeros(n,n,N);
xs = zeros(n, N);
Ps = zeros(n,n,N);

[xp(:,1), Pp(:,:,1)] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
[xf(:,1), Pf(:,:,1)] = nonLinKFupdate(xp(:,1), Pp(:,:,1), Y(:,1), S(:,1), h, R, sigmaPoints, type);

    for i = 2:1:N
        [xp(:,i), Pp(:,:,i)] = nonLinKFprediction(xf(:,i-1), Pf(:,:,i-1), f, T, Q,sigmaPoints, type);
        [xf(:,i), Pf(:,:,i)] = nonLinKFupdate(xp(:,i), Pp(:,:,i), Y(:,i), S(:,i), h, R, sigmaPoints, type);
    end
    
    xs(:,N) = xf(:,N); 
    Ps(:,:,N) = Pf(:,:,N);
    for i = N-1:-1:1
        [xs(:,i), Ps(:,:,i)] = nonLinRTSSupdate(xs(:,i+1), Ps(:,:,i+1), xf(:,i), Pf(:,:,i), xp(:,i+1), Pp(:,:,i+1), f, T, sigmaPoints, type);
    end
    %function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, f, T, sigmaPoints, type)
    
end

function [x, P] = nonLinKFprediction(x, P, f, T, Q, sigmaPoints, type)
    % Your code here! Basically the same as in HA3, but with a function 
    % handle for the sigmaPoints generating function instead of a
    % direct call.
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

function [x, P] = nonLinKFupdate(x, P, y, s, h, R, sigmaPoints, type)
    % Your code here! Basically the same as in HA3, but with a function 
    % handle for the sigmaPoints generating function instead of a
    % direct call.
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

function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, ...
                                     f, T, sigmaPoints, type)
    % Your code here! Copy from previous task!
    switch type
    case 'EKF'
        [~, A] = f(xf_k, T);
        
        G_k = Pf_k*A'*Pp_kplus1^-1;
        xs = xf_k + G_k*(xs_kplus1-f(xf_k, T));
        Ps = Pf_k - G_k*(Pp_kplus1 - Ps_kplus1)*G_k';
    otherwise
        [SP,W] = sigmaPoints(xf_k, Pf_k, type);
        P_k_kplus1 = zeros(length(xf_k));
        for i = 1:length(W)
            P_k_kplus1 = P_k_kplus1 + (SP(:,i) - xf_k)*(f(SP(:,i),T) - xp_kplus1)'*W(i);
        end
        
        G_k = P_k_kplus1*Pp_kplus1^-1;
        xs = xf_k + G_k*(xs_kplus1 - xp_kplus1);
        Ps = Pf_k - G_k*(Pp_kplus1 - Ps_kplus1)*G_k';
    end
end