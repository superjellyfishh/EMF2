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
% T = 100 for the last point
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
    X(1:end,1) = delta_p;
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
histogram(tstat_ar1,30,'FaceColor',[0 0.4470 0.7410]) % NOT centered around zero!
title("t-stats of beta in (2.2)")

% Dickey Fuller critical values
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
X(1:end,1) = diff(log(DAT(:,2)));
X(1:end,2) = DAT(2:end,2);
LM_US = fitlm(X(:,1),X(:,2))
beta = LM_US.Coefficients{2,1};
tstat_ar1_US = (LM_US.Coefficients{2,1})/(LM_US.Coefficients{2,2}) 
% Can't reject : the prices of stocks in the US are difference stationary
% (DS)

% For UK Stock Market
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) =  diff(log(DAT(:,4)));
X(1:end,2) = DAT(2:end,4);
LM_UK = fitlm(X(:,1),X(:,2))
beta = LM_UK.Coefficients{2,1};
tstat_ar1_UK = (LM_UK.Coefficients{2,1})/(LM_UK.Coefficients{2,2}) 
% Can't reject : the prices of stocks in the US are difference stationary
% (DS)

% We need to convert the dividend yield into dividend payments D:
Div_Pay_US = (DAT(:,2) .* (DAT(:,3))/100) / 12;
Div_Pay_UK = (DAT(:,4) .* (DAT(:,5))/100) / 12;

% For US Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = Div_Pay_US(2:end);
X(1:end,2) = Div_Pay_US(1:end-1);
LM_D_US = fitlm(X(:,1),X(:,2))
beta = LM_D_US.Coefficients{2,1};
tstat_ar1_D_US = (LM_D_US.Coefficients{2,1} - 1)/(LM_D_US.Coefficients{2,2})
% We reject H0 -> not stationary

% For UK Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = Div_Pay_UK(2:end);
X(1:end,2) = Div_Pay_UK(1:end-1);
LM_D_UK = fitlm(X(:,1),X(:,2))
beta = LM_D_UK.Coefficients{2,1};
tstat_ar1_D_UK = (LM_D_UK.Coefficients{2,1} - 1)/(LM_D_UK.Coefficients{2,2})


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
    % Es timate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = p(1:end-1);
    X(1:end,2) = p(2:end);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_96(i) = (LM.Coefficients{2,1} - 1)/(LM.Coefficients{2,2});
    
    disp(i)
end

histogram(tstat_ar1_96, 'FaceColor',[0 0.4470 0.7410])
title('Distribution of t-stats of beta of (2.5)')
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
power10_80 = sum(tstat_ar1_80 < c10)/length(tstat_ar1_80)

% At 5%:
power5_80 = sum(tstat_ar1_80 < c5)/length(tstat_ar1_80)

% At 1%:
power1_80 = sum(tstat_ar1_80 < c1)/length(tstat_ar1_80)

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
T = 360;
N = 10000;
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
    LM_z = fitlm(Xz(:,1),Xz(:,2));
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

% Sames values!

%% Exercise 3b
% Estimate the regression for the US and UK markets:
% p(t) = a + b*d(t) + z(t)
T = length(DAT(:,2));
% But first we need to convert the dividend yield into dividend payments D:
Div_Pay_US = (DAT(:,2) .* (DAT(:,3))/100) / 12;
Div_Pay_UK = (DAT(:,4) .* (DAT(:,5))/100) / 12;

% For the US Market:
LM_US = fitlm(log(DAT(:,2)),log(Div_Pay_US))

% Dickey-Fuller test for z in the US:
z = LM_US.Residuals{:,1};
Xz = zeros(T-1,2);
Xz(:,1) = z(1:end-1);
Xz(1:end,2) = z(2:end);
LM_z = fitlm(Xz(:,1),Xz(:,2));
tstat_ar1_zUS = (LM_z.Coefficients{2,1} - 1)/(LM_z.Coefficients{2,2});
% T-stat = -1.7476

% For the UK Market:
LM_UK = fitlm(log(DAT(:,4)),log(Div_Pay_UK))

% Dickey-Fuller test for z in the US:
z = LM_UK.Residuals{:,1};
Xz = zeros(T-1,2);
Xz(:,1) = z(1:end-1);
Xz(1:end,2) = z(2:end);
LM_z = fitlm(Xz(:,1),Xz(:,2));
tstat_ar1_zUK = (LM_z.Coefficients{2,1} - 1)/(LM_z.Coefficients{2,2});
% T-stat = -1.836

% Comment your results and conclude on the cointegration/non- 
% cointegration between the stock price and dividend processes. 
% Looking at the values of the parameter estimates (a and b), give
% your conclusion about the DDM.

%% 3c. Error-correcting model for UK

delta_p = diff(log(DAT(:,4)));
delta_pt1 = delta_p(2:end);
delta_dt1 = diff(log(Div_Pay_UK));
zt1 = z(2:end);
Y = [delta_pt1 delta_dt1(2:end) zt1(2:end)];

LM_err_UK = fitlm(delta_p(2:end),Y) % Doesn't work................
