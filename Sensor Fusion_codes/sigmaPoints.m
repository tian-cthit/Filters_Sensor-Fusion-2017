function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] matrix with sigma points
%   W           [1 x 2n+1] vector with sigma point weights 
%
switch type        
    case 'UKF'
        % your code
        n = length(x);
        w_0 = 1 - n/3;
        sqrt_P = sqrtm(P);
        
        SP = zeros(n, 2*n + 1);
        SP(:,1) = x;
        for i = 2:1:n+1
                SP(:,i) = x + sqrt(3)*sqrt_P(:,i-1);
        end
        for i = n+2:1:2*n+1
                SP(:,i) = x - sqrt(3)*sqrt_P(:,i-n-1);
        end
        
        
        W = zeros(1, 2*n + 1);
        W(1) = w_0;        
        W(2:2*n+1) = (1-w_0)/(2*n);        
         
        
    case 'CKF'
        n = length(x);
        sqrt_P = sqrtm(P);
        w_0 = 0;
        
        SP = zeros(n, 2*n);
        
        for i = 1:1:n
            SP(:,i) = x + sqrt(n/(1-w_0))*sqrt_P(:,i);
        end
        for i = n+1:1:2*n
            SP(:,i) = x - sqrt(n/(1-w_0))*sqrt_P(:,i-n);
        end
        W = zeros(1, 2*n);
        for i = 1:1:2*n
            W(i) = (1-w_0)/(2*n);
        end
        
    otherwise
        error('Incorrect type of sigma point')
end
end