function [xfp, Pfp, Xp, Wp, j] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.

n = size(x_0,1);
m = size(Y,1);
K = size(Y,2);

xfp = zeros(n,K);
Pfp = zeros(n,n,K);
Xp = zeros(n,N,K);
Wp = zeros(N,K);

Xp_0 = mvnrnd(x_0,P_0,N)';
W_0 = ones(1,N)/N;

[X_k, W_k] = pfFilterStep(Xp_0, W_0, Y(:,1), proc_f, proc_Q, meas_h, meas_R);  
% plotFunc(1, X_k, W_k, xf, Pf, bResample, sigma, ax);

hold on;
if bResample == true
   [X_k, W_k, j] = resampl(X_k, W_k); 
end

Xp(:,:,1) = X_k;
Wp(:,1) = W_k';

xfp(:,1) = sum(X_k.*W_k,2);
Pfp(:,:,1) = (X_k - xfp(:,1).*ones(n,N))*((X_k - xfp(:,1).*ones(n,N)).*W_k)';

for k = 2:1:K
    [X_k, W_k] = pfFilterStep(Xp(:,:,k-1), Wp(:,k-1)', Y(:,k), proc_f, proc_Q, meas_h, meas_R);
%     plotFunc(k, X_k, W_k, xf, Pf, bResample, sigma, ax);
    
    if bResample == true
    [X_k, W_k, j] = resampl(X_k, W_k); 
    end
    
    Xp(:,:,k) = X_k;
    Wp(:,k) = W_k'; 
    
    xfp(:,k) = sum(X_k.*W_k,2);
    Pfp(:,:,k) = (X_k - xfp(:,k).*ones(n,N))*((X_k - xfp(:,k).*ones(n,N)).*W_k)';
end
    

end

function [Xr, Wr, j] = resampl(X, W)
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


function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, y_k, proc_f, proc_Q, meas_h, meas_R)

n = size(X_kmin1,1);
N = size(X_kmin1,2);
X_k = zeros(n,N);
Y_k_pre = zeros(length(y_k),1);

    X_k = proc_f(X_kmin1(:,:)) + mvnrnd(zeros(1,n),proc_Q,N)';
    
    Y_k_pre = meas_h(X_k);    %m x N    


Y = mvnpdf(y_k',Y_k_pre',meas_R);  % Y:Nx1
W_k_0 = W_kmin1 .* Y';   
W_k_sum = sum(W_k_0);
W_k = W_k_0/W_k_sum;
end