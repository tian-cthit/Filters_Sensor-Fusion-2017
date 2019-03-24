clear
clc
mu_x = [10;0];
Sigma_x = [0.2 0;0 8];
N = 5000; %samples
A = [1 1;1 -1];
b = 0;
%% 2.1a
R = mvnrnd(mu_x,Sigma_x,N);         %sampling
Sample = cell(1,N);
for i = 1:1:N
    Sample{i} = A*(R(i,:)');        %samples go through function
end
ESum = zeros(length(Sample{1}),1);
CSum = zeros(length(Sample{1}),length(Sample{1}));
for i = 1:1:N
    ESum = ESum + Sample{i};        %mean of samples
end
E = ESum/N;
for i = 1:1:N
    CSum = CSum + (Sample{i} - E)*(Sample{i} - E)';         % covariance of samples
end
C = CSum/N;
mu_y = E;
Sigma_y = C;
y = zeros(N,2);
for i = 1:1:N
y(i,:) = Sample{i};     
end
y_s = y;

%%
[mu_y2, Sigma_y2] = affineGaussianTransform(mu_x, Sigma_x, A, b);   % calculate mean and covariance 
%% 2.1a plot

figure
for i = 1:1:5000
    h1 = plot(Sample{i}(1),Sample{i}(2),'.k');  %plot samples
    hold on;
end

hold on;
[ a, p1] = sigmaEllipse2D( mu_y, Sigma_y);      %plot 3-sigma ellipse
[ b, p2] = sigmaEllipse2D( mu_y2, Sigma_y2);
h2 = plot(mu_y(1),mu_y(2),'b+');                %plot mean
h3 = plot(mu_y2(1),mu_y2(2),'r+');
legend([h1 p1 p2 h2 h3],'+Samples','3-Sigma region by approxGaussianTransform','3-Sigma region by affineGaussianTransform','mean by approxGaussianTransform','mean by affineGaussianTransform');

%%
R = mvnrnd(mu_x,Sigma_x,N);   % sampling
Sample = cell(1,N);
for i = 1:1:N
    Sample{i} = [norm(R(i,:));atan2(R(i,2),R(i,1))];        %samples go through function
end

ESum = zeros(length(Sample{1}),1);
CSum = zeros(length(Sample{1}),length(Sample{1}));
for i = 1:1:N
    ESum = ESum + Sample{i};       
end
E = ESum/N;     % mean
for i = 1:1:N
    CSum = CSum + (Sample{i} - E)*(Sample{i} - E)';
end
C = CSum/N;     % covariance
mu_y = E;
Sigma_y = C;
y = zeros(N,2); 
for i = 1:1:N
y(i,:) = Sample{i}; 
end
y_s = y;
%% 2.1b plot
figure

for i = 1:1:N
    h1 = plot(Sample{i}(1),Sample{i}(2),'.k');  %plot samples
    hold on;
end
axis([8, 12, -1, 1.4]);
hold on;
[ a, p1] = sigmaEllipse2D( mu_y, Sigma_y);      %plot 3-sigma ellipse

h2 = plot(mu_y(1),mu_y(2),'b+');                %plot mean

 legend([h1 p1 h2],'+Samples','3-Sigma region by approxGaussianTransform','mean by approxGaussianTransform');
%%




