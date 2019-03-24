clear
clc
%% True track
% Sampling period
T  = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(200:400) = -pi/201/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
% Simulate
X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
% Set turn rate
X(5,i) = omega(i);
end
%%
N = 600;
S = [280*ones(1,N); -140*ones(1,N)];
sigma_r = 15;
sigma_phi = 2*pi/180;
R = [sigma_r^2 0;0 sigma_phi^2];
h = @rangeBearingMeasurements;
f = @coordinatedTurnMotion;

Y = genNonLinearMeasurementSequence(X, S, h, R);

%% a
x_cart = zeros(1,N);
y_cart = zeros(1,N);
for i = 1:1:N
     [x_cart(i), y_cart(i)] = pol2cart(Y(2,i),Y(1,i));
     x_cart(i) = x_cart(i) + 280;
     y_cart(i) = y_cart(i) - 140;
end

figure(1) %state 1
hold on;
measurements1 = plot(1:N,x_cart(:),'r');
states1 = plot(X(1,:),'b');
xlabel('T');
ylabel('position');
legend([measurements1,states1],'measurements','true states');

figure(2) %state 2
hold on;
measurements2 = plot(1:N,y_cart(:),'r');
states2 = plot(X(2,:),'b');
xlabel('T');
ylabel('position');
legend([measurements2,states2],'measurements','true states');
%% b
sigma_v = [1,0.01,5,10,20,0];
sigma_omega = [pi/180, 0.01*pi/180, 4*pi/180, 8*pi/180, 16*pi/180, 0];

P_0 = 0.001*diag([1 1 1 1 1]);
%% 
flag = 4;
switch flag
    case 1
        sigma_vv = sigma_v(1);
        sigma_oo = sigma_omega(1);
    case 2
        sigma_vv = sigma_v(2);
        sigma_oo = sigma_omega(2);
    case 3
        sigma_vv = sigma_v(1);
        sigma_oo = sigma_omega(2);
    case 4
        sigma_vv = sigma_v(2);
        sigma_oo = sigma_omega(1);
    case 5
        sigma_vv = 0.005;
        sigma_oo = 0.8*pi/180;
end
Q = diag([0; 0; sigma_vv^2; 0; sigma_oo^2]);
x_0 = [0 0 20 0 omega(1)]';

[xf_UKF,Pf_UKF,xp_UKF,Pp_UKF] = nonLinearKalmanFilter(Y,x_0,P_0,f,T,Q,S,h,R,'UKF');
%%
figure(1) %state 1
hold on;
measurements1 = plot(1:N,x_cart(:),'r');
states1 = plot(X(1,:),'k');
pred1 = plot(xp_UKF(1,:),'b');
est1 = plot(xf_UKF(1,:),'c');
xlabel('T');
ylabel('position');
legend([measurements1,states1,pred1,est1],'measurements','true states','prediction','estimation');

figure(2) %state 2
hold on;
measurements2 = plot(1:N,y_cart(:),'r');
states2 = plot(X(2,:),'k');
pred2 = plot(xp_UKF(2,:),'b');
est2 = plot(xf_UKF(2,:),'c');
xlabel('T');
ylabel('position');
legend([measurements2,states2,pred2,est2],'measurements','true states','prediction','estimation');
%%
figure(3) %state 3
hold on;
states3 = plot(X(3,:),'k');
pred3 = plot(xp_UKF(3,:),'b');
est3 = plot(xf_UKF(3,:),'r');
xlabel('T');
ylabel('position');
legend([states3,pred3,est3],'true states','prediction','estimation');

figure(4) %state 4
hold on;
states4 = plot(X(4,:),'k');
pred4 = plot(xp_UKF(4,:),'b');
est4 = plot(xf_UKF(4,:),'r');
xlabel('T');
ylabel('position');
legend([states4,pred4,est4],'true states','prediction','estimation');

figure(5) %state 5
hold on;
states5 = plot([1:N],X(5,2:N+1),'k');
pred5 = plot(xp_UKF(5,:),'b');
est5 = plot(xf_UKF(5,:),'r');
xlabel('T');
ylabel('position');
legend([states5,pred5,est5],'true states','prediction','estimation');
%%
figure(6)
hold on;
a = plot(X(1,:),X(2,:),'k');
b = plot(x_cart,y_cart,'--b');
c = plot(xp_UKF(1,:),xp_UKF(2,:),'r');
xlabel('x');
ylabel('y');
legend([a,b,c],'true states','measurement','estimation');