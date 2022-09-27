%% data prepare
analysis.path='D:\MATLAB\R2021a\workspace\BCD';
analysis.hspicepath = '"C:\synopsys\Hspice_B-2007.03\BIN\hspice.exe"';
cd(analysis.path);
% parameters set
[analysis.model.mean,analysis.model.sigma]=initial_format();
mu=analysis.model.mean;
sigma=analysis.model.sigma;
% train sample 
load('M_train.mat','M_X_train','M_Y_train');

%% normalization
HowManySamples=length(M_X_train);

M_X_train=(M_X_train-repmat(mu,HowManySamples,1))./repmat(sigma,HowManySamples,1);
trmin=min(M_Y_train); trmax=max(M_Y_train);
M_Y_train=(M_Y_train-repmat(trmin,HowManySamples,1))./(trmax-trmin);

H=M_X_train';
Z_samples=M_Y_train;
%% calculated weights

%% polynomial lasso
% tmp=ones(1,1000);
H2=H';
D=x2fx(H2,'linear');
[coef,coef0]=polylasso(H,Z_samples);
Y_train_lasso=D*coef+coef0;
Z_samples=M_Y_train-Y_train_lasso;

% sparsity selection
[~,indcol]=find(coef>0);
H(indcol)=[];
%% RBF interpolation
r=1000;
K=@(x,y)exp(-1/r*norm(x-y)^2);

tic
GramMatrix2=exp(-1/r*(repmat(diag(H'*H)',HowManySamples,1) ...
    -2*(H')*H ...
    + repmat(diag(H'*H)',HowManySamples,1)'));
toc

Weights=pinv(GramMatrix2)*Z_samples;

%% evaluating the approximation

%% normalization
load('M_test.mat','M_X_test','M_Y_test');
HowManyTestSamples=length(M_X_test);
M_X_test=(M_X_test-repmat(mu,HowManyTestSamples,1))./repmat(sigma,HowManyTestSamples,1);
M_Y_test=(M_Y_test-repmat(trmin,HowManyTestSamples,1))./(trmax-trmin);

%% lasso polynomial evaluation
D_t=x2fx(M_X_test,'linear');
Z_appro_part1=D_t*coef+coef0;

%% RBF interp evaluation
BigEval=M_X_test';
Z_appro_part2=Weights'*exp(-1/r*(repmat(diag(H'*H),1,size(BigEval,2))-2*H'*BigEval+repmat(diag(BigEval'*BigEval)',HowManySamples,1)));
% Z_approximation02=reshape(Z_approximation,length(t_evaluate),[]);        

Z_appro_part2=Z_appro_part2';
% Z_appro=Z_appro_part1+Z_appro_part2;

% Z_approximation=Z_approximation';
%% compare
Er=immse(M_Y_test,Z_approximation);
plot(M_Y_test,'*');
hold on
plot(Z_approximation,'*');

figure
hist(M_Y_test,16);
figure
hist(Z_approximation,16)