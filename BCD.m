%% data prepare
clc;clear;
analysis.path='D:\MATLAB\R2021a\workspace\BCD';
analysis.hspicepath = '"C:\synopsys\Hspice_B-2007.03\BIN\hspice.exe"';
cd(analysis.path);
% parameters set
[analysis.model.mean,analysis.model.sigma]=initial_format();
mu=analysis.model.mean;
sigma=analysis.model.sigma;
% train sample 

% M_X_train01=M_X_train(:,3);
% M_X_train02=M_X_train(:,3);

%% normalization
load('M_train.mat','M_X_train','M_Y_train');
HowManySamples=length(M_X_train);

M_X_train=(M_X_train-repmat(mu,HowManySamples,1))./repmat(sigma,HowManySamples,1);
trmin=min(M_Y_train); trmax=max(M_Y_train);
M_Y_train=(M_Y_train-repmat(trmin,HowManySamples,1))./(trmax-trmin);

% H=M_X_train';
% Z_samples=M_Y_train;

Add_model=zeros(length(M_Y_train),1);
Res=(M_Y_train-Add_model);
FRES=sum(abs(Res))/length(Res);
[row,col]=size(M_X_train); unit=3;
iteration=0; M_X_train_tmp=[];
% Weights_old=[]; 
% Z_appro_part_old=zeros(length(M_Y_train),0);

while iteration<100
   
    k=0; Weights=[]; Z_appro=[];
%     Z_appro_old=Z_appro;

    

for i=1:col/unit
    k=k+1;
    M_X_train_tmp=M_X_train(:,unit*(i-1)+1:unit*i);
%     HowManySamples=length(M_X_train_tmp);
    H=M_X_train_tmp';
    Z_samples=M_Y_train;
    %% polynomial lasso
% % tmp=ones(1,1000);
    H2=H';
    D=x2fx(H2,'linear');
    [coef,coef0]=polylasso(D,Z_samples);
    Y_train_lasso=D*coef+coef0;
    Z_samples=M_Y_train-Y_train_lasso;

%     % sparsity selection
%     [~,indcol]=find(coef>0);
%     H(indcol)=[];    
%% RBF part    
%     Z_samples=M_Y_train;
    r=1000;
    
    % training 
    K=@(x,y)exp(-1/r*norm(x-y)^2);
    % calculate weights
    tic
    GramMatrix2=exp(-1/r*(repmat(diag(H'*H)',HowManySamples,1) ...
    -2*(H')*H ...
    + repmat(diag(H'*H)',HowManySamples,1)'));
    toc

    Weights_tmp=pinv(GramMatrix2)*Z_samples;
    Weights=[Weights,Weights_tmp];
    
    %% resdual check by test set 
    % normalization
    load('M_test.mat','M_X_test','M_Y_test');
    HowManyTestSamples=length(M_X_test);
    M_X_test=(M_X_test-repmat(mu,HowManyTestSamples,1))./repmat(sigma,HowManyTestSamples,1);
    M_Y_test=(M_Y_test-repmat(trmin,HowManyTestSamples,1))./(trmax-trmin);

    M_X_test=M_X_test(:,unit*(i-1)+1:unit*i);
    % calc res
    %% lasso polynomial evaluation
    D_t=x2fx(M_X_test,'linear');
    Z_appro_part1=D_t*coef+coef0;
    
    %% RBF interp evaluation
    BigEval=M_X_test';
    Z_appro_part2=Weights_tmp'*exp(-1/r*(repmat(diag(H'*H),1,size(BigEval,2))-2*H'*BigEval+repmat(diag(BigEval'*BigEval)',HowManySamples,1)));
    % Z_approximation02=reshape(Z_approximation,length(t_evaluate),[]);        

    Z_appro_part2=Z_appro_part2';
    Z_appro_tmp=Z_appro_part1+Z_appro_part2;

%     % Z_approximation=Z_approximation';
%     BigEval=M_X_train_tmp';
%     Z_appro_part2=Weights_tmp'*exp(-1/r*(repmat(diag(H'*H),1,size(BigEval,2))-2*H'*BigEval+repmat(diag(BigEval'*BigEval)',HowManySamples,1)));
% % Z_approximation02=reshape(Z_approximation,length(t_evaluate),[]);        
%     Z_appro_part2=Z_appro_part2';
    Z_appro=[Z_appro,Z_appro_tmp];
    Res=Res-(Z_appro(k)-Z_appro(k));
    FRES=sum(abs(Res));
end
    iteration=iteration+1;
    B=cumsum(Z_appro,2);
    Add_model=B(:,end);
    Res=M_Y_test-Add_model;
    FRES=sum(abs(Res))/length(Res);

end


%% Residual updating