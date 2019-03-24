clear
clc
%% a
N = 100;
x_0 = [0;0;14;0;0];
P_0 = diag([10^2 10^2 2^2 (pi/180)^2 (5*pi/180)^2]);
S = 100*ones(2,N);

f = @coordinatedTurnMotion;
h = @rangeBearingMeasurements;

T = 1;
%% 
flag = 2;
if flag == 1
    Q = diag([0 0 1 0 pi/180]);
    R = diag([25 pi/180]);
else
    Q = diag([0 0 1 0 pi/180]);
    R = diag([1 10*pi/180]);
end

X = genNonLinearStateSequence(x_0, P_0, f, T, Q, N);
Y = genNonLinearMeasurementSequence(X, S, h, R);
[xf_EKF,Pf_EKF,xp_EKF,Pp_EKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'EKF');
[xf_UKF,Pf_UKF,xp_UKF,Pp_UKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'UKF');
[xf_CKF,Pf_CKF,xp_CKF,Pp_CKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'CKF');

x_cart = zeros(1,N);
y_cart = zeros(1,N);
for i = 1:1:N
     [x_cart(i), y_cart(i)] = pol2cart(Y(2,i),Y(1,i));
     x_cart(i) = x_cart(i) + 100;
     y_cart(i) = y_cart(i) + 100;
end
%%  state space plot
% figure(1);
% hold on;
% plot(X(1,:),X(2,:),'b');
% plot(x_cart(:),y_cart(:),'r');
%%
for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    v1(i) = Pf_EKF(1,1,i);    %transfer P in into an array
    v2(i) = Pf_EKF(2,2,i); 
end
figure(2)   % state 1 EKF
hold on;
a2ekf = plot([0:N],X(1,:),'K');
b2ekf = plot([1:N],x_cart(:),'r+');
c2ekf = plot([1:N], [xf_EKF(1,:) ], 'b');
d2ekf = plot([1:N], [xf_EKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ekf = plot([1:N], [xf_EKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ekf,b2ekf,c2ekf,d2ekf,e2ekf],'true states','measurements','EKF estimation','+3-sigma region','-3-sigma region');

figure(3)   % state 1 uKF
hold on;
a2ukf = plot([0:N],X(1,:),'K');
b2ukf = plot([1:N],x_cart(:),'r+');
c2ukf = plot([1:N], [xf_UKF(1,:) ], 'b');
d2ukf = plot([1:N], [xf_UKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ukf = plot([1:N], [xf_UKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ukf,b2ukf,c2ukf,d2ukf,e2ukf],'true states','measurements','UKF estimation','+3-sigma region','-3-sigma region');

figure(4)   % state 1 cKF
hold on;
a2ckf = plot([0:N],X(1,:),'K');
b2ckf = plot([1:N],x_cart(:),'r+');
c2ckf = plot([1:N], [xf_CKF(1,:) ], 'b');
d2ckf = plot([1:N], [xf_CKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ckf = plot([1:N], [xf_CKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ckf,b2ckf,c2ckf,d2ckf,e2ckf],'true states','measurements','CKF estimation','+3-sigma region','-3-sigma region');

%%
figure(5)   % state 2 EKF
hold on;
a2ekf = plot([0:N],X(2,:),'K');
b2ekf = plot([1:N],y_cart(:),'r+');
c2ekf = plot([1:N], [xf_EKF(2,:) ], 'b');
d2ekf = plot([1:N], [xf_EKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ekf = plot([1:N], [xf_EKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ekf,b2ekf,c2ekf,d2ekf,e2ekf],'true states','measurements','EKF estimation','+3-sigma region','-3-sigma region');

figure(6)   % state 1 uKF
hold on;
a2ukf = plot([0:N],X(2,:),'K');
b2ukf = plot([1:N],y_cart(:),'r+');
c2ukf = plot([1:N], [xf_UKF(2,:) ], 'b');
d2ukf = plot([1:N], [xf_UKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ukf = plot([1:N], [xf_UKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ukf,b2ukf,c2ukf,d2ukf,e2ukf],'true states','measurements','UKF estimation','+3-sigma region','-3-sigma region');

figure(7)   % state 1 cKF
hold on;
a2ckf = plot([0:N],X(2,:),'K');
b2ckf = plot([1:N],y_cart(:),'r+');
c2ckf = plot([1:N], [xf_CKF(2,:) ], 'b');
d2ckf = plot([1:N], [xf_CKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ckf = plot([1:N], [xf_CKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ckf,b2ckf,c2ckf,d2ckf,e2ckf],'true states','measurements','CKF estimation','+3-sigma region','-3-sigma region');
%%  compute the mean square error
error_ekf1 = xf_EKF(1:2,:) - X(1:2,2:N+1);
mean_cov_ekf = zeros(2,1);
error_ukf1 = xf_UKF(1:2,:) - X(1:2,2:N+1);
mean_cov_ukf = zeros(2,1);
error_ckf1 = xf_CKF(1:2,:) - X(1:2,2:N+1);
mean_cov_ckf = zeros(2,1);
for i = 1:1:N
    mean_cov_ekf = mean_cov_ekf + error_ekf1(:,i) * error_ekf1(:,i)';
    mean_cov_ukf = mean_cov_ukf + error_ukf1(:,i) * error_ukf1(:,i)';
    mean_cov_ckf = mean_cov_ckf + error_ckf1(:,i) * error_ckf1(:,i)';
end
mean_cov_ekf = mean_cov_ekf/N;
mean_cov_ukf = mean_cov_ukf/N;
mean_cov_ckf = mean_cov_ckf/N;
%% prediction plot
for i = 1:N
    %error_cov(i) = (X(i) - X_F(i))*(X(i) - X_F(i))';  % error covariance
    v1(i) = Pp_EKF(1,1,i);    %transfer P in into an array
    v2(i) = Pp_EKF(2,2,i); 
end
figure(12)   % state 1 EKF
hold on;
a2ekf = plot([0:N],X(1,:),'K');
b2ekf = plot([1:N],x_cart(:),'r+');
c2ekf = plot([1:N], [xp_EKF(1,:) ], 'b');
d2ekf = plot([1:N], [xp_EKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ekf = plot([1:N], [xp_EKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ekf,b2ekf,c2ekf,d2ekf,e2ekf],'true states','measurements','EKF prediction','+3-sigma region','-3-sigma region');

figure(13)   % state 1 uKF
hold on;
a2ukf = plot([0:N],X(1,:),'K');
b2ukf = plot([1:N],x_cart(:),'r+');
c2ukf = plot([1:N], [xp_UKF(1,:) ], 'b');
d2ukf = plot([1:N], [xp_UKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ukf = plot([1:N], [xp_UKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ukf,b2ukf,c2ukf,d2ukf,e2ukf],'true states','measurements','UKF prediction','+3-sigma region','-3-sigma region');

figure(14)   % state 1 cKF
hold on;
a2ckf = plot([0:N],X(1,:),'K');
b2ckf = plot([1:N],x_cart(:),'r+');
c2ckf = plot([1:N], [xp_CKF(1,:) ], 'b');
d2ckf = plot([1:N], [xp_CKF(1,:) ] + 3*sqrt(v1(:)'), '--b');
e2ckf = plot([1:N], [xp_CKF(1,:) ] - 3*sqrt(v1(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ckf,b2ckf,c2ckf,d2ckf,e2ckf],'true states','measurements','CKF prediction','+3-sigma region','-3-sigma region');

%%
figure(15)   % state 2 EKF
hold on;
a2ekf = plot([0:N],X(2,:),'K');
b2ekf = plot([1:N],y_cart(:),'r+');
c2ekf = plot([1:N], [xp_EKF(2,:) ], 'b');
d2ekf = plot([1:N], [xp_EKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ekf = plot([1:N], [xp_EKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ekf,b2ekf,c2ekf,d2ekf,e2ekf],'true states','measurements','EKF prediction','+3-sigma region','-3-sigma region');

figure(16)   % state 1 uKF
hold on;
a2ukf = plot([0:N],X(2,:),'K');
b2ukf = plot([1:N],y_cart(:),'r+');
c2ukf = plot([1:N], [xp_UKF(2,:) ], 'b');
d2ukf = plot([1:N], [xp_UKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ukf = plot([1:N], [xp_UKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ukf,b2ukf,c2ukf,d2ukf,e2ukf],'true states','measurements','UKF prediction','+3-sigma region','-3-sigma region');

figure(17)   % state 1 cKF
hold on;
a2ckf = plot([0:N],X(2,:),'K');
b2ckf = plot([1:N],y_cart(:),'r+');
c2ckf = plot([1:N], [xp_CKF(2,:) ], 'b');
d2ckf = plot([1:N], [xp_CKF(2,:) ] + 3*sqrt(v2(:)'), '--b');
e2ckf = plot([1:N], [xp_CKF(2,:) ] - 3*sqrt(v2(:)'), '--b');
xlabel('T');
ylabel('position');
legend([a2ckf,b2ckf,c2ckf,d2ckf,e2ckf],'true states','measurements','CKF prediction','+3-sigma region','-3-sigma region');
%% d
MC = 2;
error_est_ekf = cell(1,1,MC);
error_est_ukf = cell(1,1,MC);
error_est_ckf = cell(1,1,MC);

error_pred_ekf = cell(1,1,MC);
error_pred_ukf = cell(1,1,MC);
error_pred_ckf = cell(1,1,MC);
for i = 1:MC
X = genNonLinearStateSequence(x_0, P_0, f, T, Q, N);
Y = genNonLinearMeasurementSequence(X, S, h, R);
[xf_EKF,Pf_EKF,xp_EKF,Pp_EKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'EKF');
[xf_UKF,Pf_UKF,xp_UKF,Pp_UKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'UKF');
[xf_CKF,Pf_CKF,xp_CKF,Pp_CKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'CKF');

error_est_ekf{i} = xf_EKF - X(:,2:N+1);
error_est_ukf{i} = xf_UKF - X(:,2:N+1);
error_est_ckf{i} = xf_CKF - X(:,2:N+1);

error_pred_ekf{i} = xp_EKF - X(:,2:N+1);
error_pred_ukf{i} = xp_UKF - X(:,2:N+1);
error_pred_ckf{i} = xp_CKF - X(:,2:N+1);    
end
%%
error_est_ekf_aver = zeros(5,N);
error_est_ukf_aver = zeros(5,N);
error_est_ckf_aver = zeros(5,N);

error_pred_ekf_aver = zeros(5,N);
error_pred_ukf_aver = zeros(5,N);
error_pred_ckf_aver = zeros(5,N);
for i = 1:MC
    error_est_ekf_aver = error_est_ekf_aver + error_est_ekf{i};
    error_est_ukf_aver = error_est_ukf_aver + error_est_ukf{i};
    error_est_ckf_aver = error_est_ckf_aver + error_est_ckf{i};
    
    error_pred_ekf_aver = error_pred_ekf_aver + error_pred_ekf{i};
    error_pred_ukf_aver = error_pred_ukf_aver + error_pred_ukf{i};
    error_pred_ckf_aver = error_pred_ckf_aver + error_pred_ckf{i};
end
error_est_ekf_aver = error_est_ekf_aver/MC;
error_est_ukf_aver = error_est_ukf_aver/MC;
error_est_ckf_aver = error_est_ckf_aver/MC;

error_pred_ekf_aver = error_pred_ekf_aver/MC;
error_pred_ukf_aver = error_pred_ukf_aver/MC;
error_pred_ckf_aver = error_pred_ckf_aver/MC;
%%  est EKF
figure(8)
histogram(error_est_ekf_aver(1,:),'Normalization','probability','BinWidth',20); 
figure(9)
histogram(error_est_ekf_aver(2,:),'Normalization','probability','BinWidth',10); 
figure(10)
histogram(error_est_ekf_aver(3,:),'Normalization','probability','BinWidth',2); 
figure(11)
histogram(error_est_ekf_aver(4,:),'Normalization','probability','BinWidth',10); 
figure(12)
histogram(error_est_ekf_aver(5,:),'Normalization','probability','BinWidth',0.3); 
%% UKF
figure(8)
histogram(error_est_ukf_aver(1,:),'Normalization','probability','BinWidth',3); 
figure(9)
histogram(error_est_ukf_aver(2,:),'Normalization','probability','BinWidth',3); 
figure(10)
histogram(error_est_ukf_aver(3,:),'Normalization','probability','BinWidth',1); 
figure(11)
histogram(error_est_ukf_aver(4,:),'Normalization','probability','BinWidth',0.2); 
figure(12)
histogram(error_est_ukf_aver(5,:),'Normalization','probability','BinWidth',0.1); 
%% CKF
figure(8)
histogram(error_est_ckf_aver(1,:),'Normalization','probability','BinWidth',10); 
figure(9)
histogram(error_est_ckf_aver(2,:),'Normalization','probability','BinWidth',8); 
figure(10)
histogram(error_est_ckf_aver(3,:),'Normalization','probability','BinWidth',1.5); 
figure(11)
histogram(error_est_ckf_aver(4,:),'Normalization','probability','BinWidth',10); 
figure(12)
histogram(error_est_ckf_aver(5,:),'Normalization','probability','BinWidth',0.5); 
%%  pred EKF
figure(8)
histogram(error_pred_ekf_aver(1,:),'Normalization','probability','BinWidth',20); 
figure(9)
histogram(error_pred_ekf_aver(2,:),'Normalization','probability','BinWidth',10); 
figure(10)
histogram(error_pred_ekf_aver(3,:),'Normalization','probability','BinWidth',2); 
figure(11)
histogram(error_pred_ekf_aver(4,:),'Normalization','probability','BinWidth',10); 
figure(12)
histogram(error_pred_ekf_aver(5,:),'Normalization','probability','BinWidth',0.3); 
%% UKF
figure(8)
histogram(error_pred_ukf_aver(1,:),'Normalization','probability','BinWidth',3); 
figure(9)
histogram(error_pred_ukf_aver(2,:),'Normalization','probability','BinWidth',3); 
figure(10)
histogram(error_pred_ukf_aver(3,:),'Normalization','probability','BinWidth',1); 
figure(11)
histogram(error_pred_ukf_aver(4,:),'Normalization','probability','BinWidth',0.2); 
figure(12)
histogram(error_pred_ukf_aver(5,:),'Normalization','probability','BinWidth',0.1); 
%% CKF
figure(8)
histogram(error_pred_ckf_aver(1,:),'Normalization','probability','BinWidth',10); 
figure(9)
histogram(error_pred_ckf_aver(2,:),'Normalization','probability','BinWidth',8); 
figure(10)
histogram(error_pred_ckf_aver(3,:),'Normalization','probability','BinWidth',1.5); 
figure(11)
histogram(error_pred_ckf_aver(4,:),'Normalization','probability','BinWidth',10); 
figure(12)
histogram(error_pred_ckf_aver(5,:),'Normalization','probability','BinWidth',0.5);

