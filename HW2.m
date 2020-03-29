clear, clc
%% First we need to import the data
D = importdata('DATA_HW2.xlsx');
DAT = D.data.Sheet1;
% Converting the dates in matlab readable dates...
DAT(:,1) = DAT(:,1) + 693960;

%% Exercise 2a : T = 360

% Simulate N time series of T error terms
N = 10000;
T = 360;
beta = zeros(N,1);
tstat_ar1 = zeros(N,1);

for i = 1:N
    err = randn(T,1);
    % err = normrnd(0,1,[1,T]);
    % Compute a time series of stock price p(t) = P(t-1) + e(t)
    p = zeros(T,1);
    
    for j = 2:T
        p(j) = p(j-1) + err(j);
    end
    delta_p = diff(p);
    alpha = 0;
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = p(1:end-1);
    X(1:end,2) = p(2:end);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1(i) = (LM.Coefficients{2,1}-1)/(LM.Coefficients{2,2});
    
%     sum1 = 0;
%     sum2 = 0;
%     for j = 2:T
%         sum1 = sum1 + p(j-1)*err(j);
%         sum2 = sum2 + p(j-1).^2;
%     end
% 
%     tstat_ar1(i) = ((T^-1)*(sum1))/(((T^-2))*(sum2));

% 	temp = ar((p),1);
% 	beta(i) = temp.Report.Parameters.ParVector;
% 	tstat_ar1(i) = ((temp.Report.Parameters.ParVector) - 1) / ...
%        (sqrt(temp.Report.Parameters.FreeParCovariance));
   
%     AR_p = arima('Constant',NaN,'ARLags',1,'Distribution','Gaussian');
%     [AR_p,EstParamCov] = estimate(AR_p,p,'Display','off');
%     tstat_ar1(i) = (AR_p.AR{1,1} - 1)/sqrt(EstParamCov(2,2));
    disp(i)
    % X = [ones(N-(testLags+1),1),yLags((testLags+2):end,2),deltaYLags((testLags+2):end,2:end)];
    % [a(i),b(i),c(i),d(i)] = adftest(p,'lags',1,'model','ARD');
end
%tstat_ar1 = (beta - 1)/std(beta);

% Histogram of the t-stats
hist(tstat_ar1,30) % NOT centered around zero!
% We also observe that there are a lot of cases when
% we rejet H0 based just on noise.

%% Dickey Fuller critical values
t_sorted = sort(tstat_ar1);

% 10%
c10 = t_sorted(N/10)

% 5%
c5 = t_sorted(N/20)

% 1%
c1 = t_sorted(N/100)

% Same values!

%% Exercise 2b
% For US Stock Market:
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = DAT(1:end-1,2);
X(1:end,2) = DAT(2:end,2);
LM_US = fitlm(X(:,1),X(:,2))
beta = LM_US.Coefficients{2,1};
tstat_ar1_US = (LM_US.Coefficients{2,1}-1)/(LM_US.Coefficients{2,2}) 
% Can't reject : the prices of stocks in the US are difference stationary
% (DS)

% For UK Stock Market
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = DAT(1:end-1,4);
X(1:end,2) = DAT(2:end,4);
LM_UK = fitlm(X(:,1),X(:,2))
beta = LM_UK.Coefficients{2,1};
tstat_ar1_UK = (LM_UK.Coefficients{2,1}-1)/(LM_UK.Coefficients{2,2}) 
% Can't reject : the prices of stocks in the US are difference stationary
% (DS)

% For US Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = DAT(1:end-1,3);
X(1:end,2) = DAT(2:end,3);
LM_D_US = fitlm(X(:,1),X(:,2))
beta = LM_D_US.Coefficients{2,1};
tstat_ar1_D_US = (LM_D_US.Coefficients{2,1}-1)/(LM_D_US.Coefficients{2,2})
% We can reject H0 at 10% confidence

% For UK Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = DAT(1:end-1,5);
X(1:end,2) = DAT(2:end,5);
LM_D_UK = fitlm(X(:,1),X(:,2))
beta = LM_D_UK.Coefficients{2,1};
tstat_ar1_D_UK = (LM_D_UK.Coefficients{2,1}-1)/(LM_D_UK.Coefficients{2,2})
% We can reject H0 at 10% confidence (was close)

% -> Dividend processes could be stationary in levels

%% 2c
N = 10000;
T = 360;
beta = zeros(N,1);
tstat_ar1_96 = zeros(N,1);

for i = 1:N
    % Compute a time series of stock price p(t) = 0.96 * P(t-1) + e(t)
    err = randn(T,1);
    p = zeros(T,1);
    for j = 2:T
        p(j) = 0.96 * p(j-1) + err(j);
    end
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = p(1:end-1);
    X(1:end,2) = p(2:end);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_96(i) = (LM.Coefficients{2,1}-1)/(LM.Coefficients{2,2});
    
    disp(i)
    
    % std_err = std(p(1:end-1))/sqrt(length(p(1:end-1)));
    % tstat_ar1_96(i) = (beta - 1)/std_err;
end

hist(tstat_ar1_96)
% Probability of rejecting H0 given that H1 is true:
% At 10%:
power10 = sum(tstat_ar1_96 < c10)/length(tstat_ar1_96)

% At 5%:
power5 = sum(tstat_ar1_96 < c5)/length(tstat_ar1_96)

% At 1%:
power1 = sum(tstat_ar1_96 < c1)/length(tstat_ar1_96)

% Values of power are pretty low considering that H1 is true..

% Plot together the cumulative distribution function of the t-stats under 
% the null (computed in 2a)and the cumulative distribution function of 
% the t-stats under the AR(1) process given above. Comment.
[counts1, bins1] = histcounts(tstat_ar1);
[counts2, bins2] = histcounts(tstat_ar1_96);
cdf1 = cumsum(counts1);
cdf2 = cumsum(counts2);
plot(cdf1), hold on
plot(cdf2)
title("CDF of t-stats")
legend("AR(1) of 2a","AR(1) of 2c")
% The t-stats are close to the distribution, when they should be very 
% shifted to the negative, so the test is .questionnable.

%% 2c: Same case but with T = 100
N = 10000;
T = 100;
beta = zeros(N,1);
tstat_ar1_T100 = zeros(N,1);

for i = 1:N
    % Compute a time series of stock price p(t) = 0.96 * P(t-1) + e(t)
    err = randn(T,1);
    p = zeros(T,1);
    for j = 2:T
        p(j) = 0.96 * p(j-1) + err(j);
    end
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = p(1:end-1);
    X(1:end,2) = p(2:end);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_T100(i) = (LM.Coefficients{2,1}-1)/(LM.Coefficients{2,2});
    
    disp(i)
    
    % std_err = std(p(1:end-1))/sqrt(length(p(1:end-1)));
    % tstat_ar1_96(i) = (beta - 1)/std_err;
end

hist(tstat_ar1_T100)
% Probability of rejecting H0 given that H1 is true:
% At 10%:
power10_T100 = sum(tstat_ar1_96 < c10)/length(tstat_ar1_T100)

% At 5%:
power5_T100 = sum(tstat_ar1_96 < c5)/length(tstat_ar1_T100)

% At 1%:
power1_T100 = sum(tstat_ar1_96 < c1)/length(tstat_ar1_T100)

% Values of power are pretty low considering that H1 is true..

% Plot together the cumulative distribution function of the t-stats under 
% the null (computed in 2a)and the cumulative distribution function of 
% the t-stats under the AR(1) process given above. Comment.
[counts1, bins1] = histcounts(tstat_ar1);
[counts2, bins2] = histcounts(tstat_ar1_T100);
cdf1 = cumsum(counts1);
cdf2 = cumsum(counts2);
plot(cdf1), hold on
plot(cdf2)
title("CDF of t-stats")
legend("AR(1) of 2a","AR(1) of 2c when T = 100")


%% 2c: Same case but with AR(1) : p(t) = 0.8*p(t-1) + err
N = 10000;
T = 360;
beta = zeros(N,1);
tstat_ar1_80 = zeros(N,1);

for i = 1:N
    % Compute a time series of stock price p(t) = 0.80 * P(t-1) + e(t)
    err = randn(T,1);
    p = zeros(T,1);
    for j = 2:T
        p(j) = 0.80 * p(j-1) + err(j);
    end
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = p(1:end-1);
    X(1:end,2) = p(2:end);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_80(i) = (LM.Coefficients{2,1}-1)/(LM.Coefficients{2,2});
    
    disp(i)
end

hist(tstat_ar1_80)
% Probability of rejecting H0 given that H1 is true:
% At 10%:
power10_3 = sum(tstat_ar1_80 < c10)/length(tstat_ar1_80)

% At 5%:
power5_3 = sum(tstat_ar1_80 < c5)/length(tstat_ar1_80)

% At 1%:
power1_3 = sum(tstat_ar1_80 < c1)/length(tstat_ar1_80)

% Values of power are much better

[counts1, bins1] = histcounts(tstat_ar1);
[counts2, bins2] = histcounts(tstat_ar1_80);
cdf1 = cumsum(counts1);
cdf2 = cumsum(counts2);
plot(cdf1), hold on
plot(cdf2)
title("CDF of t-stats")
legend("AR(1) of 2a","AR(1) of 2c (case when phi(1) = 0.8)")

% Seems to be better when the AR(1) process is further away from 1

%% Exercise 3a
T = 100;
N = 1000;
beta_3a = zeros(N,1);
tstat_3a = zeros(N,1);
tstat_ar1_z = zeros(N,1);
z = zeros(N,1);

for i = 1:N
    err = randn(T,2);

    p = zeros(T,1);
    d = zeros(T,1);
    for j = 2:T
        p(j) = p(j-1) + err(j,1);
        d(j) = d(j-1) + err(j,2);
    end
    X = zeros(T,2);
    X(1:end,1) = p(1:end);
    X(1:end,2) = d(1:end);
    LM_3a3 = fitlm(X(:,1),X(:,2));
    z = LM_3a3.Residuals{:,1};
    
    Xz = zeros(T-1,2);
    Xz(:,1) = z(1:end-1);
    Xz(1:end,2) = z(2:end);
    LM_z = fitlm(X(:,1),X(:,2));
    tstat_ar1_z(i) = (LM_z.Coefficients{2,1} - 1)/(LM_z.Coefficients{2,2});
    disp(i)
end
hist(tstat_ar1_z)
t_sorted = sort(tstat_ar1_z);

% 10%
c10z = t_sorted(N/10)

% 5%
c5z = t_sorted(N/20)

% 1%
c1z = t_sorted(N/100)

