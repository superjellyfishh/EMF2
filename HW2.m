clear, clc
%% First we need to import the data
D = importdata('DATA_HW2.xlsx');
DAT = D.data.Sheet1;
% Converting the dates in matlab readable dates...
DAT(:,1) = DAT(:,1) + 693960;

%% Exercise 2a: T = 360

% Simulate N time series of T error terms
N = 10000;
T = 360;
beta = zeros(N,1);
tstat_ar1 = zeros(N,1);

for i = 1:N
    err = randn(T,1);
    % Compute a time series of stock price p(t) = P(t-1) + e(t)
    p = zeros(T,1);
    
    for j = 2:T
        p(j) = p(j-1) + err(j);
    end
    
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = diff(p);
    X(1:end,2) = p(1:end-1);
    LM = fitlm(X(:,1), X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1(i) = (LM.Coefficients{2,1})/(LM.Coefficients{2,2});
    disp(i)
end

% Histogram of the t-stats
histogram(tstat_ar1,30,'FaceColor',[0 0.4470 0.7410])
% NOT centered around zero!
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

%% Exercise 2a : T = 100

% Simulate N time series of T error terms
N = 10000;
T = 100;
beta = zeros(N,1);
tstat_ar1t100 = zeros(N,1);

for i = 1:N
    err = randn(T,1);
    % Compute a time series of stock price p(t) = P(t-1) + e(t)
    p = zeros(T,1);
    
    for j = 2:T
        p(j) = p(j-1) + err(j);
    end
    
    % Estimate the AR(1) model, Compute the t-stat for beta
    X = zeros(T-1,2);
    X(1:end,1) = diff(p);
    X(1:end,2) = p(1:end-1);
    LM = fitlm(X(:,1), X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1t100(i) = (LM.Coefficients{2,1})/(LM.Coefficients{2,2});
    disp(i)
end

histogram(tstat_ar1t100)
% Dickey Fuller critical values
t100_sorted = sort(tstat_ar1t100);

% 10%
c10t100 = t100_sorted(N/10)

% 5%
c5t100 = t100_sorted(N/20)

% 1%
c1t100 = t100_sorted(N/100)

%% Exercise 2b
% For US Stock Market:
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = diff(log(DAT(:,2)));
X(1:end,2) = log(DAT(1:end-1,2));
LM_US = fitlm(X(:,1),X(:,2))
beta = LM_US.Coefficients{2,1};
tstat_ar1_US = (LM_US.Coefficients{2,1} - 1)/(LM_US.Coefficients{2,2}) 
% Can't reject : the returns of stocks in the US are difference stationary
% (DS)

% For UK Stock Market:
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = diff(log(DAT(:,4)));
X(1:end,2) = log(DAT(1:end-1,4));
LM_UK = fitlm(X(:,1),X(:,2))
beta = LM_UK.Coefficients{2,1};
tstat_ar1_UK = (LM_UK.Coefficients{2,1})/(LM_UK.Coefficients{2,2}) 
% Can't reject : the returns of stocks in the US are difference stationary
% (DS)

% We need to convert the dividend yield into dividend payments D:
Div_Pay_US = (DAT(:,2) .* (DAT(:,3))/100) / 12;
Div_Pay_UK = (DAT(:,4) .* (DAT(:,5))/100) / 12;

% For US Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = diff(log(Div_Pay_US));
X(1:end,2) = log(Div_Pay_US(1:end-1));
LM_D_US = fitlm(X(:,1),X(:,2))
beta = LM_D_US.Coefficients{2,1};
tstat_ar1_D_US = (LM_D_US.Coefficients{2,1})/(LM_D_US.Coefficients{2,2})
% We reject H0 -> not stationary

% For UK Dividend Process
T = size(DAT,1);
X = zeros(T-1,2);
X(1:end,1) = diff(log(Div_Pay_UK));
X(1:end,2) = log(Div_Pay_UK(1:end-1));
LM_D_UK = fitlm(X(:,1),X(:,2))
beta = LM_D_UK.Coefficients{2,1};
tstat_ar1_D_UK = (LM_D_UK.Coefficients{2,1})/(LM_D_UK.Coefficients{2,2})

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
    X(1:end,1) = diff(p(1:end));
    X(1:end,2) = p(1:end-1);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_96(i) = (LM.Coefficients{2,1})/(LM.Coefficients{2,2});
    
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
% shifted to the negative, so the test is questionnable.

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
    X(1:end,1) = p(2:end);
    X(1:end,2) = p(1:end-1);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_T100(i) = (LM.Coefficients{2,1} - 1)/(LM.Coefficients{2,2});
    
    disp(i)
end

histogram(tstat_ar1_T100)
% Probability of rejecting H0 given that H1 is true:
% At 10%:
power10_T100 = sum(tstat_ar1_96 < c10t100)/length(tstat_ar1_T100)

% At 5%:
power5_T100 = sum(tstat_ar1_96 < c5t100)/length(tstat_ar1_T100)

% At 1%:
power1_T100 = sum(tstat_ar1_96 < c1t100)/length(tstat_ar1_T100)

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
    X(1:end,1) = diff(p(1:end));
    X(1:end,2) = p(1:end - 1);
    LM = fitlm(X(:,1),X(:,2));
    beta(i) = LM.Coefficients{2,1};
    tstat_ar1_80(i) = (LM.Coefficients{2,1})/(LM.Coefficients{2,2});
    
    disp(i)
end

histogram(tstat_ar1_80)
% Probability of rejecting H0 given that H1 is true:
% At 10%:
power10_80 = sum(tstat_ar1_80 < c10)/length(tstat_ar1_80)

% At 5%:
power5_80 = sum(tstat_ar1_80 < c5)/length(tstat_ar1_80)

% At 1%:
power1_80 = sum(tstat_ar1_80 < c1)/length(tstat_ar1_80)

% Values of power are much better

[counts1, bins1] = histcounts(tstat_ar1);
[counts2c, bins2] = histcounts(tstat_ar1_80);
cdf1 = cumsum(counts1);
cdf2c = cumsum(counts2c);
plot(cdf1), hold on
plot(cdf2c)
title("CDF of t-stats")
legend("AR(1) of 2a","AR(1) of 2c (case when phi(1) = 0.8)")

% Seems to be better when the AR(1) process is further away from 1
histogram(tstat_ar1, 'FaceColor',[0 0.4470 0.7410])
hold on
histogram(tstat_ar1_96, 'FaceColor','r')
histogram(tstat_ar1_T100)
histogram(tstat_ar1_80)
title("CDF of t-stats")
legend("AR(1) of 2a","AR(1) phi(1) = 0.96)", ...
    "AR(1) T = 100)", ...
    "AR(1) phi(1) = 0.8)")
xlim([-9 4])

%% Exercise 3a
T = 360;
% T = 100 for comparing the results in the remark
N = 10000;
beta_3a = zeros(N,1);
tstat_3a = zeros(N,1);
tstat_ar1_z = zeros(N,1);
% z = zeros(N,1);

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
    Xz(:,1) = diff(z);
    Xz(1:end,2) = z(1:end-1);
    LM_z = fitlm(Xz(:,1),Xz(:,2));
    tstat_ar1_z(i) = (LM_z.Coefficients{2,1})/(LM_z.Coefficients{2,2});
    disp(i)
end
histogram(tstat_ar1_z,'FaceColor',[0.4660 0.6740 0.1880])
title('t-stats of beta of (3.2)')
t_sorted = sort(tstat_ar1_z);

% 10%
c10z = t_sorted(N/10)

% 5%
c5z = t_sorted(N/20)

% 1%
c1z = t_sorted(N/100)

% Sames values when T = 100

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
Xz(:,1) = diff(z);
Xz(:,2) = z(1:end - 1);
LM_z = fitlm(Xz(:,1),Xz(:,2));
tstat_ar1_zUS = (LM_z.Coefficients{2,1})/(LM_z.Coefficients{2,2})
% T-stat = -1.7476

% For the UK Market:
LM_UK = fitlm(log(DAT(:,4)),log(Div_Pay_UK))

% Dickey-Fuller test for z in the US:
z = LM_UK.Residuals{:,1};
Xz = zeros(T-1,2);
Xz(:,1) = diff(z);
Xz(1:end,2) = z(1:end - 1);
LM_z = fitlm(Xz(:,1),Xz(:,2));
tstat_ar1_zUK = (LM_z.Coefficients{2,1})/(LM_z.Coefficients{2,2})
% T-stat = -1.836

% Comment your results and conclude on the cointegration/non- 
% cointegration between the stock price and dividend processes. 
% Looking at the values of the parameter estimates (a and b), give
% your conclusion about the DDM.
plot(cumsum(diff(log(DAT(:,2)))))
hold on
plot(cumsum(diff(log(Div_Pay_US))))
title("Stock Prices & Dividend Payments in the US")
legend("Stock Prices","Dividend Payments")

%% 3c. Error-correcting model for UK

UKp = log(DAT(:,4));
UKd = log((DAT(:,4).*(DAT(:,5))/100)/12);
UKzTemp = fitlm(UKd, UKp);
UKz = UKzTemp.Residuals{:,1};

delta_UKp = diff(UKp);
delta_UKd = diff(UKd);
X = [delta_UKp(1:(end-1),:), delta_UKd(1:(end-1),:), UKz(2:(end-1),:)];
y_1 = delta_UKp(2:end,:);
y_2 = delta_UKd(2:end,:);

pUK = fitlm(X, y_1);
dUK = fitlm(X, y_2);

VECMp = pUK.Coefficients;
VECMp.Properties.RowNames = {'\mu_1','\varphi_{11}','\varphi_{12}','\gamma_1'};
VECMd = dUK.Coefficients;
VECMd.Properties.RowNames = {'\mu_2','\varphi_{21}','\varphi_{22}','\gamma_2'};
VECM = [VECMp;VECMd]
