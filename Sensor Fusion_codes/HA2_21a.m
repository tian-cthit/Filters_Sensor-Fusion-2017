clear
clc
%%
N = 20;
n = 1;
% Motion Parameters
A = 1;
Q = 1.5;
% Measurement Parameters
H = 1;
R = 2.5;
% Prior
xPrior = 0;
PPrior = 2.6;
%% a)
%generate states & measurements
X = genLinearStateSequence(xPrior, PPrior, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);

%% plots
plot(0:1:N,X);
hold on
plot(1:N,Y);
xlabel('Samples');
ylabel('Value');
legend('States','Measurements');
%% b)
% filtered data, where :
% X_F: filtered state
% P: filtered cov
% Xp: prediction states
% Pp = prediction cov
[X_F, P, Xp, Pp] = kalmanFilter_with_prediction(Y, xPrior, PPrior, A, Q, H, R);

%% plots
figure(2)
hold on;
plot(X, '*k'); % real states
plot(Y, '*r'); % measurements

% 3-sigma region
plot([0:N], [xPrior X_F], 'b');
plot([0:N], [xPrior X_F] + 3*sqrt([PPrior P(:)']), '--b');
plot([0:N], [xPrior X_F] - 3*sqrt([PPrior P(:)']), '--b');
title('Reference solution');
xlabel('k');
ylabel('x');
legend('real states','measurements', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');

%% error covariance
error_cov = zeros(1,N);
for i = 1:N
    error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    p(i) = P(:,:,i);    %transfer P in into an array
end
e = error_cov - p;  %error between error cov and uncertainty
figure(3)
plot(1:N,e);
xlabel('k');
ylabel('e');
legend('Error between error covariance and uncertainty');
%% c
figure(4)
grid on;
hold on;
plot(0:N,X, '*r'); % real states
plot([1:N], X_F, '*k');

% posterior density
mu = zeros(1,N);
sigma2 = zeros(1,N);
for i = 1:3:N
    %[mu(i), sigma2(i)] = posteriorGaussian(X_F(i), p(i), Y(i), R);
    k = 1;
    for j = X_F(i)-10:0.01:X_F(i)+10
        posterior(k) = normpdf(j,X_F(i),P(i));
        k = k+1;
    end
    plot(posterior+i,X_F(i)-10:0.01:X_F(i)+10,'b')
end
xlabel('k');
ylabel('e');
legend('Real states','Filtered states','Posterior density');
%% convariance & uncertainty
% error_cov_c = zeros(1,N);
% e_c = mu - p;  %error between error cov and uncertainty
% figure(5)
% plot(1:N,e);
% xlabel('k');
% ylabel('e');
% legend('Error between error covariance and uncertainty');
%% d
figure(8)
hold on;
meas = plot(Y, '*b'); % measurements
% posterior at k:  p(x_k|y1:k)
for i = 1:1:N
    k = 1;
    for j = X_F(i)-10:0.01:X_F(i)+10
        posterior(k) = normpdf(j,X_F(i),p(i));
        k = k+1;
    end
    p1 = plot(posterior+i,X_F(i)-10:0.01:X_F(i)+10,'b');   
end

% posterior at k-1: p(x_k-1|y1:k-1)
for i = 1:1:N
    k = 1;
    for j = X_F(i)-10:0.01:X_F(i)+10
        posterior(k) = normpdf(j,X_F(i),p(i));
        k = k+1;
    end
    p2 = plot(posterior+i+1,X_F(i)-10:0.01:X_F(i)+10,'k')   ;
end

% p(x_k|y1:k-1)
for i = 1:1:N
    k = 1;
    for j = Xp(i)-10:0.01:Xp(i)+10
        prediction(k) = normpdf(j,Xp(i),Pp(i));
        k = k+1;
    end
    p3 = plot(prediction+i-1,Xp(i)-10:0.01:Xp(i)+10,'r')   ;
end
xlabel('Probability+k');
ylabel('state value');
legend([meas p1,p2,p3], 'measurements','posterior at k','posterior at k-1','prediction at k-1');
%% e
N_e = 500;
%% 
%generate states & measurements
X_e = genLinearStateSequence(xPrior, PPrior, A, Q, N_e);
Y_e = genLinearMeasurementSequence(X_e, H, R);
% filtered data
[X_F_e, P_e] = kalmanFilter(Y_e, xPrior, PPrior, A, Q, H, R);

e_e = X_e(2:N_e+1) - X_F_e;
%% 
figure(6)
histogram(e_e,round(max(e_e)-min(e_e)),'Normalization','probability')
hold on;
plot(-20:0.1:20,normpdf(-20:0.1:20,0,P_e(N)))
xlabel('Mean Error');
ylabel('Probability');
legend('Mean Error','N(x;0,P(N))');
%% f
N = 20;
xPrior_f = 10;
[X_F_f, P_f] = kalmanFilter(Y, xPrior_f, PPrior, A, Q, H, R);

figure(7)
hold on;
plot(X, '.k'); % real states
plot([0:N], [xPrior X_F], 'r');
plot([0:N], [xPrior X_F_f], 'b');
legend('Real States','Correct Mean','Incorrect Mean');
xlabel('Samples');
ylabel('Value');