clear
clc
%% 
x1_0 = [500;0];
p1_0 = [50^2 0;0 1];
x2_0 = [100;-50];
p2_0 = [200 180;180 200];
x3_0 = [500;0];
p3_0 = [1 0;0 50^2];

sigma_r = 1;
sigma_p = pi/180;
sigma_measurement = [sigma_r 0;0 pi/180];
s = [0;0];
%% a

N = 8000;
x1 = zeros(2,N);
x2 = zeros(2,N);
x3 = zeros(2,N);
y1 = zeros(2,N);
y2 = zeros(2,N);
y3 = zeros(2,N);
h = @rangeBearingMeasurements;
for i = 1:1:N
    x1(:,i) = mvnrnd(x1_0,p1_0)';
    y1(:,i) = h(x1(:,i),s) + mvnrnd([0;0],sigma_measurement)';
    x2(:,i) = mvnrnd(x2_0,p2_0)';
    y2(:,i) = h(x2(:,i),s) + mvnrnd([0;0],sigma_measurement)';
    x3(:,i) = mvnrnd(x3_0,p3_0)';
    y3(:,i) = h(x3(:,i),s) + mvnrnd([0;0],sigma_measurement)';
end
y1_mean = [sum(y1(1,:));sum(y1(2,:))]/N;
y2_mean = [sum(y2(1,:));sum(y2(2,:))]/N;
y3_mean = [sum(y3(1,:));sum(y3(2,:))]/N;
y1_sub = (y1 - y1_mean);
y2_sub = (y2 - y2_mean);
y3_sub = (y3 - y3_mean);
y1_cov = zeros(2,2);
y2_cov = zeros(2,2);
y3_cov = zeros(2,2);
for i = 1:1:N
    y1_cov = y1_cov + y1_sub(:,i)*(y1_sub(:,i))';
end
    y1_cov = y1_cov/N;
for i = 1:1:N
    y2_cov = y2_cov + y2_sub(:,i)*(y2_sub(:,i))';
end
    y2_cov = y2_cov/N;
for i = 1:1:N
    y3_cov = y3_cov + y3_sub(:,i)*(y3_sub(:,i))';
end
    y3_cov = y3_cov/N;
%%
figure(1);
hold on;
a1 = plot(y1(1,:),y1(2,:),'b.');
b1 = plot(y1_mean(1),y1_mean(2),'r+');
[ aaa, p1] = sigmaEllipse2D( y1_mean, y1_cov);
xlabel('y(1)');
legend([a1,b1,p1],'samples','mean','3-sigma ellipse');

figure(2);
hold on;
a2 = plot(y2(1,:),y2(2,:),'b.');
b2 = plot(y2_mean(1),y2_mean(2),'r+');
[ aaa2, p2] = sigmaEllipse2D( y2_mean, y2_cov);
xlabel('y(1)');
legend([a2,b2,p2],'samples','mean','3-sigma ellipse');

figure(3);
hold on;
a3 = plot(y3(1,:),y3(2,:),'b.');
b3 = plot(y3_mean(1),y3_mean(2),'r+');
[ aaa3, p3] = sigmaEllipse2D( y3_mean, y3_cov);
xlabel('y(1)');
legend([a3,b3,p3],'samples','mean','3-sigma ellipse');
%% b
h = @rangeBearingMeasurements;
R = sigma_measurement;

% EKF
[x_ekf1, A] = h(x1_0,s);
P_ekf1 = A*p1_0*A' + R;
[x_ekf2, A] = h(x2_0,s);
P_ekf2 = A*p2_0*A' + R;
[x_ekf3, A] = h(x3_0,s);
P_ekf3 = A*p3_0*A' + R;
%%
% UKF
n = length(x1_0);
%1
[SP_ukf1,W_ukf1] = sigmaPoints(x1_0, p1_0, 'UKF');        
x_ukf1 = zeros(n,1);
for i = 1:1:2*n+1
    x_ukf1 = x_ukf1 + h(SP_ukf1(:,i),s)*W_ukf1(i);
end
P_ukf1 = zeros(n,n);
for i = 1:1:2*n+1
    P_ukf1 = P_ukf1 + (h(SP_ukf1(:,i),s) - x_ukf1)*(h(SP_ukf1(:,i),s) - x_ukf1)'*W_ukf1(i);
end
P_ukf1 = P_ukf1 + R;  

%2        
[SP_ukf2,W_ukf2] = sigmaPoints(x2_0, p2_0, 'UKF');        
x_ukf2 = zeros(n,1);
for i = 1:1:2*n+1
    x_ukf2 = x_ukf2 + h(SP_ukf2(:,i),s)*W_ukf2(i);
end
P_ukf2 = zeros(n,n);
for i = 1:1:2*n+1
    P_ukf2 = P_ukf2 + (h(SP_ukf2(:,i),s) - x_ukf2)*(h(SP_ukf2(:,i),s) - x_ukf2)'*W_ukf2(i);
end
P_ukf2 = P_ukf2 + R; 

% 3
[SP_ukf3,W_ukf3] = sigmaPoints(x3_0, p3_0, 'UKF');        
x_ukf3 = zeros(n,1);
for i = 1:1:2*n+1
    x_ukf3 = x_ukf3 + h(SP_ukf3(:,i),s)*W_ukf3(i);
end
P_ukf3 = zeros(n,n);
for i = 1:1:2*n+1
    P_ukf3 = P_ukf3 + (h(SP_ukf3(:,i),s) - x_ukf3)*(h(SP_ukf3(:,i),s) - x_ukf3)'*W_ukf3(i);
end
P_ukf3 = P_ukf3 + R; 
%% ckf
n = length(x1_0);
%1
[SP_ckf1,W_ckf1] = sigmaPoints(x1_0, p1_0, 'CKF');        
x_ckf1 = zeros(n,1);
for i = 1:1:2*n
    x_ckf1 = x_ckf1 + h(SP_ckf1(:,i),s)*W_ckf1(i);
end
P_ckf1 = zeros(n,n);
for i = 1:1:2*n
    P_ckf1 = P_ckf1 + (h(SP_ckf1(:,i),s) - x_ckf1)*(h(SP_ckf1(:,i),s) - x_ckf1)'*W_ckf1(i);
end
P_ckf1 = P_ckf1 + R;  

%2        
[SP_ckf2,W_ckf2] = sigmaPoints(x2_0, p2_0, 'CKF');        
x_ckf2 = zeros(n,1);
for i = 1:1:2*n
    x_ckf2 = x_ckf2 + h(SP_ckf2(:,i),s)*W_ckf2(i);
end
P_ckf2 = zeros(n,n);
for i = 1:1:2*n
    P_ckf2 = P_ckf2 + (h(SP_ckf2(:,i),s) - x_ckf2)*(h(SP_ckf2(:,i),s) - x_ckf2)'*W_ckf2(i);
end
P_ckf2 = P_ckf2 + R; 

% 3
[SP_ckf3,W_ckf3] = sigmaPoints(x3_0, p3_0, 'CKF');        
x_ckf3 = zeros(n,1);
for i = 1:1:2*n
    x_ckf3 = x_ckf3 + h(SP_ckf3(:,i),s)*W_ckf3(i);
end
P_ckf3 = zeros(n,n);
for i = 1:1:2*n
    P_ckf3 = P_ckf3 + (h(SP_ckf3(:,i),s) - x_ckf3)*(h(SP_ckf3(:,i),s) - x_ckf3)'*W_ckf3(i);
end
P_ckf3 = P_ckf3 + R; 
%%
figure(4);
hold on;
ekf_mean1 = plot(x_ekf1(1),x_ekf1(2),'r.');
ukf_mean1 = plot(x_ukf1(1),x_ukf1(2),'b.');
ckf_mean1 = plot(x_ckf1(1),x_ckf1(2),'k.');
b1 = plot(y1_mean(1),y1_mean(2),'+');
[aaa1, ekf_cov1] = sigmaEllipse2D(x_ekf1, P_ekf1);
[bbb1, ukf_cov1] = sigmaEllipse2D(x_ukf1, P_ukf1);
[ccc1, ckf_cov1] = sigmaEllipse2D(x_ckf1, P_ckf1);
[ aaa, p1] = sigmaEllipse2D( y1_mean, y1_cov);
xlabel('y(1)');
ylabel('y(2)');
legend([ekf_mean1, ukf_mean1,ckf_mean1,b1, ekf_cov1, ukf_cov1, ckf_cov1,p1],'EKF mean','UKF mean','CKF mean','sapmle mean','EKF 3-sigma ellipse','UKF 3-sigma ellipse','CKF 3-sigma ellipse','sample 3-sigma ellipse');

%2
figure(5);
hold on;
ekf_mean2 = plot(x_ekf2(1),x_ekf2(2),'r.');
ukf_mean2 = plot(x_ukf2(1),x_ukf2(2),'b.');
ckf_mean2 = plot(x_ckf2(1),x_ckf2(2),'k.');
b1 = plot(y2_mean(1),y2_mean(2),'+');
[aaa2, ekf_cov2] = sigmaEllipse2D(x_ekf2, P_ekf2);
[bbb2, ukf_cov2] = sigmaEllipse2D(x_ukf2, P_ukf2);
[ccc2, ckf_cov2] = sigmaEllipse2D(x_ckf2, P_ckf2);
[ aaa, p2] = sigmaEllipse2D( y2_mean, y2_cov);
xlabel('y(1)');
ylabel('y(2)');
legend([ekf_mean2, ukf_mean2,ckf_mean2,b1, ekf_cov2, ukf_cov2, ckf_cov2,p2],'EKF mean','UKF mean','CKF mean','sapmle mean','EKF 3-sigma ellipse','UKF 3-sigma ellipse','CKF 3-sigma ellipse','sample 3-sigma ellipse');

%3
figure(6);
hold on;
ekf_mean3 = plot(x_ekf3(1),x_ekf3(2),'r.');
ukf_mean3 = plot(x_ukf3(1),x_ukf3(2),'b.');
ckf_mean3 = plot(x_ckf3(1),x_ckf3(2),'k.');
b1 = plot(y3_mean(1),y3_mean(2),'+');
[aaa3, ekf_cov3] = sigmaEllipse2D(x_ekf3, P_ekf3);
[bbb3, ukf_cov3] = sigmaEllipse2D(x_ukf3, P_ukf3);
[ccc3, ckf_cov3] = sigmaEllipse2D(x_ckf3, P_ckf3);
[ aaa, p3] = sigmaEllipse2D( y3_mean, y3_cov);
xlabel('y(1)');
ylabel('y(2)');
legend([ekf_mean3, ukf_mean3,ckf_mean3,b1, ekf_cov3, ukf_cov3, ckf_cov3,p3],'EKF mean','UKF mean','CKF mean','sapmle mean','EKF 3-sigma ellipse','UKF 3-sigma ellipse','CKF 3-sigma ellipse','sample 3-sigma ellipse');



