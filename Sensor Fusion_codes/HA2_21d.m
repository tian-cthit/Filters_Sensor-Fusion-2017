clear
clc
%% d
N = 500;
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
%%
%generate states & measurements
X = genLinearStateSequence(xPrior, PPrior, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);

%plots
% plot(0:1:N,X);
% hold on
% plot(1:N,Y);
% xlabel('Samples');
% ylabel('Value');
% legend('States','Measurements');

% filtered data
[X_F, P] = kalmanFilter(Y, xPrior, PPrior, A, Q, H, R);

e = X(2:N+1) - X_F;
%%
figure(1)
histogram(e,round(max(e)-min(e)),'Normalization','probability')
hold on;
plot(-20:0.1:20,normpdf(-20:0.1:20,0,P(N)))
xlabel('Mean Error');
ylabel('Probability');
legend('Mean Error','N(x;0,P(N))');
%%
