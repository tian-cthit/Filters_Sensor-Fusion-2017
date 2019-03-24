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
omega(200:400) =-pi/201/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;
% Create true track
for i=2:K+1
% Simulate
X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
% Set turn-rate
X(5,i) = omega(i);
end
%%
N = size(X,2)-1;
f = @coordinatedTurnMotion;
h = @rangeBearingMeasurements;

P_0 = diag([10^2 10^2 2^2 (pi/180)^2 (5*pi/180)^2]);
S = [280*ones(1,N); -140*ones(1,N)];

sigma_vv = 0.005;
sigma_oo = 0.8*pi/180;
Q = diag([0; 0; sigma_vv^2; 0; sigma_oo^2]);

sigma_r = 15;
sigma_phi = 2*pi/180;
R = [sigma_r^2 0;0 sigma_phi^2];
%%
Y = genNonLinearMeasurementSequence(X, S, h, R);
%% for 2.2
Y(:,150) = [500;1.4]; 
for i = 150
     [x_cart(i), y_cart(i)] = pol2cart(Y(2,i),Y(1,i));
     x_cart(i) = x_cart(i) + 280;
     y_cart(i) = y_cart(i) - 140;
end
%%
sig = @sigmaPoints;
type = 'UKF';
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x0, P_0, f, T, Q, S, h, R, sig, type);

x_cart = zeros(1,N);
y_cart = zeros(1,N);
for i = 1:1:N
     [x_cart(i), y_cart(i)] = pol2cart(Y(2,i),Y(1,i));
     x_cart(i) = x_cart(i) + 280;
     y_cart(i) = y_cart(i) - 140;
end
%%
figure(1)
hold on;
sta = plot(X(1,:),X(2,:),'K');
me = plot(x_cart(:),y_cart(:),'c.');
fil = plot(xf(1,:),xf(2,:),'R');
smo = plot(xs(1,:),xs(2,:),'B');
for i = 1:5:N
    [ xy ,p ] = sigmaEllipse2D( xs(1:2,i), Ps(1:2,1:2,i));
    p = plot(xy(1,:),xy(2,:),'r');
end

xlabel('position X');
ylabel('position Y');
legend([sta,me,fil,smo,p],'true state','measurement','filtered state','smoothed state','3-sigma region');
%% 2.2

figure(2)
hold on;
sta = plot(X(1,:),X(2,:),'K');
me = plot(x_cart(:),y_cart(:),'b.');
fil = plot(xf(1,:),xf(2,:),'R');
smo = plot(xs(1,:),xs(2,:),'B');
for i = 1:5:N
    [ xy ,p ] = sigmaEllipse2D( xs(1:2,i), Ps(1:2,1:2,i));
    p = plot(xy(1,:),xy(2,:),'r');
end

xlabel('position X');
ylabel('position Y');
legend([sta,me,fil,smo,p],'true state','measurement','filtered state','smoothed state','3-sigma region');
%%

