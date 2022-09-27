 %% initialize
    rng(0,'twister');
    n = 1000;
    X = linspace(0,1,n)';
    X = [X,X.^2];
    y = 1 + X*[1;2] + sin(20*X*[1;-2])./(X(:,1)+1) + 0.2*randn(n,1);

 %% BCD Fitting   
    gprbcd = fitrgp(X,y,'KernelFunction','squaredexponential',...
        'FitMethod','exact','PredictMethod','bcd','BlockSize',200);
 
%% ds
rng default % For reproducibility
N = 1e4; % Number of samples
p = 1e3; % Number of features
X = randn(N,p);
beta = randn(p,1); % Multiplicative coefficients
beta0 = randn; % Additive term
y = beta0 + X*beta + randn(N,1); % Last term is noise
%% dss
B = lasso(X,y,"UseCovariance",false); % Warm up lasso for reliable timing data
tic
B = lasso(X,y,"UseCovariance",false);
timefalse = toc
