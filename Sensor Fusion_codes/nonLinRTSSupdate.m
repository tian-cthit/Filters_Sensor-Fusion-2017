function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

% Your code here.

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