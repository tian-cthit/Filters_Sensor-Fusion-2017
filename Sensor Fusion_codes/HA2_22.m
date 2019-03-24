clear
clc
%%
N =500;
T = 0.5;
A = [1 T;0 1];
H = [1 0];
R = 2;
% case d) 1
    Q_1 = [0 0;0 0.1];
% case d) 2
    Q_2 = [0 0;0 1];
% case d) 3
    Q_3 = [0 0;0 10];
% case a) & d) 4 
    Q = [0 0;0 1.5];
    
xPrior = [1;3];
PPrior = 4*eye(2);
%% a
%generate states & measurements
X = genLinearStateSequence(xPrior, PPrior, A, Q, N);
Y = genLinearMeasurementSequence(X, H, R);

%% plot states 1 & measurements
figure(1)
plot(X(1,:),'r');
hold on;
plot(2:N+1,Y,'b');
xlabel('t');
ylabel('Position');
legend('Real Position','Measurements')

% plot states 2
figure(2)
plot(X(2,:),'b')
xlabel('t');
ylabel('Velocity');
legend('Velocity')
%% b
% filtered data
% X_F: filtered state
% P: filtered cov
% Xp: prediction states
% Pp = prediction cov
[X_F, P, Xp, Pp] = kalmanFilter_with_prediction(Y, xPrior, PPrior, A, Q, H, R);

%  position plots
figure(3)
hold on;
plot(X(1,:),'r') % real states
plot(Y, '.k'); % measurements

% 3-sigma region
error_cov = zeros(1,N);
for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    p(i) = P(1,1,i);    %transfer P in into an array
end

plot([0:N], [xPrior(1) X_F(1,:)], 'b');
plot([0:N], [xPrior(1) X_F(1,:)] + 3*sqrt([PPrior(1) p(:)']), '--b');
plot([0:N], [xPrior(1) X_F(1,:)] - 3*sqrt([PPrior(1) p(:)']), '--b');
title('Reference solution');
xlabel('k');
ylabel('x');
legend('real states','measurements', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');
%%
% velocity
for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    v(i) = P(2,2,i);    %transfer P in into an array
end

figure(4)
hold on;
plot(X(2,:),'r') % real states
plot([0:N], [xPrior(2) X_F(2,:)], 'b');
plot([0:N], [xPrior(2) X_F(2,:)] + 3*sqrt([PPrior(2) v(:)']), '--b');
plot([0:N], [xPrior(2) X_F(2,:)] - 3*sqrt([PPrior(2) v(:)']), '--b');
title('Reference solution');
xlabel('k');
ylabel('v');
legend('real states', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');
%% c
figure(11)
hold on;
grid on;
meas = plot(Y, '*b'); % measurements
% posterior at k:  p(x_k|y1:k)
for i = 1:1:N
    k = 1;
    for j = X_F(1,i)-10:0.01:X_F(1,i)+10
        posterior(k) = normpdf(j,X_F(1,i),p(i));
        k = k+1;
    end
    p1 = plot(posterior+i,X_F(1,i)-10:0.01:X_F(1,i)+10,'b');   
end

% posterior at k-1: p(x_k-1|y1:k-1)
for i = 1:1:N
    k = 1;
    for j = X_F(1,i)-10:0.01:X_F(1,i)+10
        posterior(k) = normpdf(j,X_F(1,i),p(i));
        k = k+1;
    end
    p2 = plot(posterior+i+1,X_F(1,i)-10:0.01:X_F(1,i)+10,'k')   ;
end

% p(x_k|y1:k-1)
for i = 1:1:N
    k = 1;
    for j = Xp(1,i)-10:0.01:Xp(1,i)+10
        prediction(k) = normpdf(j,Xp(1,i),Pp(1,i));
        k = k+1;
    end
    p3 = plot(prediction+i-1,Xp(1,i)-10:0.01:Xp(1,i)+10,'r')   ;
end
xlabel('Probability+k');
ylabel('State value');
legend([meas,p1,p2,p3],'measurements','posterior at k','posterior at k-1','prediction at k-1');
%% d 1)
[X_F_1, P_1] = kalmanFilter(Y, xPrior, PPrior, A, Q_1, H, R);

for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    p_1(i) = P_1(2,2,i);    %transfer P in into an array
end

figure(5)
hold on;
plot(X(1,:),'r') % real states
plot([0:N], [xPrior(1) X_F_1(1,:)], 'b');
plot([0:N], [xPrior(1) X_F_1(1,:)] + 3*sqrt([PPrior(1) p_1(:)']), '--b');
plot([0:N], [xPrior(1) X_F_1(1,:)] - 3*sqrt([PPrior(1) p_1(:)']), '--b');
xlabel('k');
ylabel('Position');
legend('real states', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');

%anlysis
figure(6)
e_1 = X_F_1(1,:) - X(1,2:N+1);
histogram(e_1,20,'Normalization','probability');
xlabel('Error');
ylabel('Probability');

se_1 = (e_1*e_1')/N;
%% d 2)
[X_F_2, P_2] = kalmanFilter(Y, xPrior, PPrior, A, Q_2, H, R);

for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    p_2(i) = P_2(2,2,i);    %transfer P in into an array
end

figure(7)
hold on;
plot(X(1,:),'r') % real states
plot([0:N], [xPrior(1) X_F_2(1,:)], 'b');
plot([0:N], [xPrior(1) X_F_2(1,:)] + 3*sqrt([PPrior(1) p_1(:)']), '--b');
plot([0:N], [xPrior(1) X_F_2(1,:)] - 3*sqrt([PPrior(1) p_1(:)']), '--b');
xlabel('k');
ylabel('Position');
legend('real states', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');

%anlysis
figure(8)
e_2 = X_F_2(1,:) - X(1,2:N+1);
histogram(e_2,20,'Normalization','probability');
xlabel('Error');
ylabel('Probability');

se_2 = (e_2*e_2')/N;
%% d 3)

[X_F_3, P_3] = kalmanFilter(Y, xPrior, PPrior, A, Q_3, H, R);

for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    p_3(i) = P_3(2,2,i);    %transfer P in into an array
end

figure(9)
hold on;
plot(X(1,:),'r') % real states
plot([0:N], [xPrior(1) X_F_3(1,:)], 'b');
plot([0:N], [xPrior(1) X_F_3(1,:)] + 3*sqrt([PPrior(1) p_1(:)']), '--b');
plot([0:N], [xPrior(1) X_F_3(1,:)] - 3*sqrt([PPrior(1) p_1(:)']), '--b');
xlabel('k');
ylabel('Position');
legend('real states', 'state estimate', '+3-sigma level', '-3-sigma level','Location','southeast');

%anlysis
figure(10)
e_3 = X_F_3(1,:) - X(1,2:N+1);
histogram(e_3,20,'Normalization','probability');
xlabel('Error');
ylabel('Probability');

se_3 = (e_3*e_3')/N;

e = X_F(1,:) - X(1,2:N+1);
se = (e*e')/N;
