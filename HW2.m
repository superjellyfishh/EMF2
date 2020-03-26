clear, clc
%% First we need to import the data
D = importdata('DATA_HW2.xlsx');
DAT = D.data.Sheet1;
% Converting the dates in matlab readable dates...
DAT(:,1) = DAT(:,1) + 693960;

% Modification

%% Exercise 2

% Simulate a time series of T error terms
N = 10000;
err = randn(N,1);

% Compute a time series of stock price p(t) = P(t-1) + e(t)
p = zeros(N,1);
p(1) = 100; % Starting price of 100

for i = 2:N
    p(i) = p(i-1) + err(i);
end
plot(p)

% Estimate the AR(1) model, Compute the t-stat for beta

ar1mdl = arima(1,0,0);






