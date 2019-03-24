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

num_p = 8;
proc_f = @y_equals_to_x;
meas_h = @y_equals_to_x;

%%
X = genLinearStateSequence(xPrior, PPrior, A, Q, N);
%%
%%

sigma = 1;

Y = genLinearMeasurementSequence(X, H, R);
[X_F, P, Xp, Pp] = kalmanFilter_with_prediction(Y, xPrior, PPrior, A, Q, H, R);


% Xk
xf = X_F;
Pf = P;
bResample = true;
plotFunc_handle  = @(k, Xk, Xkmin1, Wk, j)(plotPostPdf(k, Xk, Wk, xf, Pf, bResample, sigma, ax));
[xfp, Pfp, Xp, Wp, j] = pfFilter(xPrior, PPrior, Y, proc_f, Q, meas_h, R, 30, true, plotFunc_handle);
[xfp2, Pfp2, Xp2, Wp2] = pfFilter(xPrior, PPrior, Y, proc_f, Q, meas_h, R, 50, false);
proc_f = @y_equals_to_x;
meas_h = @y_equals_to_x;
%%
figure(1)
hold on;
sta = plot(X,'k');
pf1 = plot(xfp,'b');
% pf2 = plot(xfp2,'m');
kf = plot(X_F,'r');
xlabel('t');
ylabel('state');
% legend([sta,pf1,pf2,kf],'true state','PF with resample','PF without resample','KF');
%% b
j2 = [1:50];
figure(3)
hold on;
for i = 2:1:N
    disp(i)
    plotPartTrajs(i, Xp2(:,:,i), Xp2(:,:,i-1), Wp2(i), j2)   
end

b = plot(X(2:N+1),'r*');

xlabel('k');
ylabel('state');
legend(b,'true state');
%% c

figure(4)
hold on;
for i = 2:1:N
    
    plotPartTrajs(i, Xp(:,:,i), Xp(:,:,i-1), Wp(i), j)   
end

b = plot(X(2:N+1),'r*');
xlabel('k');
ylabel('state');
legend(b,'true state');
