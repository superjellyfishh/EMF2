clear, clc
%% First we need to import the data
D = importdata('DATA_HW2.xlsx');
DAT = D.data.Sheet1;
% Converting the dates in matlab readable dates...
DAT(:,1) = DAT(:,1) + 693960;

%% Exercise 2 

% Simulate N time series of T error terms
N = 10000;
T = 100;
beta = zeros(N,1);
tstat_ar1 = zeros(N,1);

for i = 1:N
    err = randn(T,1);
    % err = normrnd(0,1,[1,T]);
    % Compute a time series of stock price p(t) = P(t-1) + e(t)
    p = zeros(T,1);
    % p(1) = err(1); % Starting price of 0
    
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
