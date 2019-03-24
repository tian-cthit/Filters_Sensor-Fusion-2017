function X = genNonLinearStateSequence(x_0, P_0, f, T, Q, N)
%GENLINEARSTATESEQUENCE generates an N-long sequence of states using a 
%    Gaussian prior and a linear Gaussian process model
%
%Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   f           Motion model function handle
%   T           Sampling time
%   Q           [n x n] Process noise covariance
%   N           [1 x 1] Number of states to generate
%
%Output:
%   X           [n x N] State vector sequence
%

% Your code here
X = zeros(length(x_0),N+1);
X(:,1) = mvnrnd(x_0,P_0)';
    for i = 2:1:N+1        
        X(:,i) = f(X(:,i-1),T) + mvnrnd(zeros(length(Q),1),Q)';
    end
end