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
    p = zeros(N,1);
    p(1) = 100; % Starting price of 100
    for j = 2:T
        p(j) = p(j-1) + err(j);
    end
    % Estimate the AR(1) model, Compute the t-stat for beta
    temp = ar(p,1);
    beta(i) = temp.A(2);
    tstat_ar1(i) = temp.A(2) / temp.Report.Parameters.FreeParCovariance;
end

plot(sort(beta))