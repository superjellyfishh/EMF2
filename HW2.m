clear, clc
%% First we need to import the data
D = importdata('DATA_HW2.xlsx');
DAT = D.data.Sheet1;
% Converting the dates in matlab readable dates...
DAT(:,1) = DAT(:,1) + 693960;

%% Exercise 2 

% Simulate a time series of T error terms
N = 1000;
T = 100;
beta = zeros(N,1);
tstat_ar1 = zeros(N,1);

for i = 1:N
    err = randn(T,1);
    % Compute a time series of stock price p(t) = P(t-1) + e(t)
    p = zeros(T,1);
    p(1) = 100; % Starting price of 100
    
    for j = 2:T
        p(j) = p(j-1) + err(j);
    end
    delta_p = diff(p);
    % Estimate the AR(1) model, Compute the t-stat for beta
    temp = ar((delta_p),1);
    beta(i) = temp.A(2);
    tstat_ar1(i) = temp.A(2) / temp.Report.Parameters.FreeParCovariance;
    disp(i)
end

% Histogram of the t-stats
hist(tstat_ar1) % Indeed centered around zero!
% We also observe that there are a lot of cases when
% we rejet H0 based just on noise.

%% Dickey Fuller test
sum1 = 0;
for i = 2:T
    sum1 = sum1 + (p(i) - beta(1)*p(i-1))^2;
end
sigma_est = (1/(T-1))*sum1;
sum2 = 0;
for i = 2:T
    sum2 = sum2 + (p(i-1))^2;
end
t_df = (beta(1)-1)/(sqrt((sum1/sum2)));
















