clear
clc
%% a
mu = 1;
Sigma = 2;

x_Hat = 1;
plot(-10:0.01:12,normpdf(-10:0.01:12,x_Hat,sqrt(Sigma)));   % MMSE
hold on;

j =0 ;
for i = -1:0.01:2                                           % MAP
    j = j+1;
    d(j) = dirac(i - 1);
    if d(j) == Inf
        d(j) = 0.3;
    end
    
end
plot(-1:0.01:2,d);
legend('MMSE','MAP')
%% b
w_1 = 0.5;
u_1 = 3;
Sigma_1 = 1;
w_2 = 0.5;
u_2 = -3;
Sigma_2 = 1;
plot(-10:0.01:10,0.5*normpdf(-10:0.01:10,u_2,sqrtm(Sigma_2)) + 0.5*normpdf(-10:0.01:10,u_1,sqrtm(Sigma_1)),'b'); % MMSE
hold on;
j =0 ;
for i = -10:0.01:10  % MAP
    j = j+1;
    d(j) = dirac(i - 3);
    d2(j) = dirac(i + 3);
    if d(j) == Inf
        d(j) = 0.4;
    end
    if d2(j) == Inf
        d2(j) = 0.4;
    end
    
end
plot(-10:0.01:10,0.5*(d+d2),'r');
legend('MMSE','MAP')
%% c
w_3 = 0.4;
u_3 = 1;
Sigma_3 = 2;
w_4 = 0.6;
u_4 = 2;
Sigma_4 = 1;
plot(-10:0.01:10,w_3*normpdf(-10:0.01:10,u_3,sqrtm(Sigma_3)) + w_4*normpdf(-10:0.01:10,u_4,sqrtm(Sigma_4)),'b'); % MMSE
hold on;
j =0 ;
for i = -10:0.01:10 %MAP
    j = j+1;
    d(j) = dirac(i - u_3);
    d2(j) = dirac(i -u_4);
    if d(j) == Inf
        d(j) = 0.7;
    end
    if d2(j) == Inf
        d2(j) = 0.7;
    end
    
end
plot(-10:0.01:10,0.5*(d+d2),'r');
legend('MMSE','MAP')
%%
plot()