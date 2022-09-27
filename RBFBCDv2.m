%% data prepare
clc;clear;
format long
analysis.path='D:\MATLAB\R2021a\workspace\BCD';
analysis.hspicepath = '"C:\synopsys\Hspice_B-2007.03\BIN\hspice.exe"';
cd(analysis.path);
% parameters set 
[analysis.model.mean,analysis.model.sigma]=initial_format();
mu=analysis.model.mean;
sigma=analysis.model.sigma;

% cut point set
cut_point_norm=zeros(1,length(mu));
cut_point=cut_point_norm.*sigma+mu;
ft=run_file(analysis.hspicepath,cut_point,1);

%% random block training
% block variable selection
block_size=3;
total_length=length(mu);
% shuffe=randperm(total_length);
shuffe=1:1:6;
block_mu_saver=[]; block_sigma_saver=[];
block_saver_in=[]; block_saver_out=[];
for i=1:total_length/block_size
    batch_size=1000;
    
    % block variables select
    block_var_mu=mu(shuffe(block_size*(i-1)+1:block_size*i));
    block_var_sigma=sigma((shuffe(block_size*(i-1)+1:block_size*i)));
    normsampling=-4+8*lhsdesign(batch_size,block_size);
    block_var=normsampling.*repmat(block_var_sigma,batch_size,1)+repmat(block_var_mu,batch_size,1);

    % replace
    total_var=repmat(cut_point,batch_size,1);
    total_var(:,shuffe(block_size*(i-1)+1:block_size*i))=block_var;
    
    % simulation
    block_res=run_file(analysis.hspicepath,total_var,batch_size)-ft;
    
    % save batch
    block_mu_saver=[block_mu_saver,block_var_mu];
    block_sigma_saver=[block_sigma_saver,block_var_sigma];
    block_saver_in=[block_saver_in,block_var];
    block_saver_out=[block_saver_out,block_res];
    
end
%% save data
save('block_saver','block_saver_in','block_saver_out','block_mu_saver','block_sigma_saver','shuffe','batch_size','block_size','total_length');

%% batch nomalization
load('block_saver','block_saver_in','block_saver_out','block_mu_saver','block_sigma_saver','shuffe','batch_size','block_size','total_length');
block_saver_in=(block_saver_in-repmat(block_mu_saver,batch_size,1))./repmat(block_sigma_saver,batch_size,1);
trmin=min(block_saver_out); trmax=max(block_saver_out);
trmin=min(trmin); trmax=max(trmax);
% block_saver_out=(block_saver_out-repmat(trmin,batch_size,1))./(trmax-trmin);

% ft=(ft-trmin)./(trmax-trmin);

%% training 
weights_saver=[];H_saver=[];coef_saver=[]; coef0_saver=[]; indrow_saver=[];
for i=1:total_length/block_size
    % load data
    M_X_train_tmp=block_saver_in(:,block_size*(i-1)+1:block_size*i);
%     HowManySamples=length(M_X_train_tmp);
    H=M_X_train_tmp';
    M_Y_train=block_saver_out(:,i);
    Z_samples=M_Y_train;
    %% polynomial lasso
% tmp=ones(1,1000);
    H2=H';
    D=x2fx(H2,'linear');
    [coef,coef0]=polylasso(D,Z_samples);
    Y_train_lasso=D*coef+coef0;
    Z_samples=M_Y_train-Y_train_lasso;

    % sparsity selection
    [indrow,~]=find(coef~=0);
    H=H(indrow-1,:);  
    
    %% RBF part    
%     Z_samples=block_saver_out(:,i); % optional
    r=2000;
    
    % training 
    K=@(x,y)exp(-1/r*norm(x-y)^2);
    % calculate weights
    tic
    GramMatrix2=exp(-1/r*(repmat(diag(H'*H)',batch_size,1) ...
    -2*(H')*H ...
    + repmat(diag(H'*H)',batch_size,1)'));
    toc

    weights=pinv(GramMatrix2)*Z_samples;
   
    weights_saver=[weights_saver,weights];
    H_saver=[H_saver;H];
    coef_saver=[coef_saver,coef];
    coef0_saver=[coef0_saver,coef0];
%     indrow_saver=[indrow_saver,indrow-1];
end

%% retraining
load('M_train.mat','M_X_train','M_Y_train');
% normalization
HowManySamples=length(M_X_train);

M_X_train=(M_X_train-repmat(mu,HowManySamples,1))./repmat(sigma,HowManySamples,1);
% trmin=min(M_Y_train); trmax=max(M_Y_train);
M_Y_train=M_Y_train-ft;
M_Y_train=(M_Y_train-repmat(trmin,HowManySamples,1))./(trmax-trmin);



%% Revise Biases 
Z_appro=0; H_len=0; Z_appro_part_saver=[];
for i=1:total_length/block_size
    tBlock_var=M_X_train(:,shuffe(block_size*(i-1)+1:block_size*i));
    
    %% lasso polynomial evaluation
    D_t=x2fx(tBlock_var,'linear');
    coef=coef_saver(:,i); 
    coef0=coef0_saver(:,i);
    Z_appro_part1=D_t*coef+coef0;
    

    
%     indrow=indrow_saver(:,i);
    BigEval=tBlock_var';
    % sparsity selection
    [indrow,~]=find(coef>0);
    BigEval=BigEval(indrow-1,:);
    weights_tmp=weights_saver(:,i);
    
%     H_len=H_len+indrow;
    H_tmp=H_saver(H_len+1:H_len+length(indrow),:);
    H_len=H_len+length(indrow);
    Z_appro_part2=weights_tmp'*exp(-1/r*(repmat(diag(H_tmp'*H_tmp),1,size(BigEval,2))-2*H_tmp'*BigEval+repmat(diag(BigEval'*BigEval)',batch_size,1)));
    
    Z_appro_part=Z_appro_part1+Z_appro_part2';
    Z_appro=Z_appro+Z_appro_part;
    Z_appro_part_saver=[Z_appro_part_saver,Z_appro_part];
end

%% LS
X=[ones(size(Z_appro_part_saver(:,1))) Z_appro_part_saver(:,1) Z_appro_part_saver(:,2)];
coef_LS=regress(M_Y_train,X);
%% MLP full connect

Mdl = fitrnet(Z_appro_part_saver,M_Y_train,"Standardize",true, ...
    "LayerSizes",[10]);



%% test
% load('M_test_MC.mat','M_X_test_MC','M_Y_test_MC');
% M_X_test=M_X_test_MC;
% M_Y_test=M_Y_test_MC;

% load('M_test_TMC.mat','M_X_test_TMC','M_Y_test_TMC');
% M_X_test=M_X_test_TMC;
% M_Y_test=M_Y_test_TMC;
% T_len=length(M_Y_test);

%  load('M_test.mat','M_X_test','M_Y_test');

load('M_train.mat','M_X_train','M_Y_train');
M_X_test=M_X_train;
M_Y_test=M_Y_train;
%% test normalization;

HowManyTestSamples=length(M_X_test);
M_X_test=(M_X_test-repmat(mu,HowManyTestSamples,1))./repmat(sigma,HowManyTestSamples,1);
M_Y_test=M_Y_test-ft;
M_Y_test=(M_Y_test-repmat(trmin,HowManyTestSamples,1))./(trmax-trmin);


%% predict
Z_appro=0; H_len=0; Z_appro_part_saver=[];
for i=1:total_length/block_size
    tBlock_var=M_X_test(:,shuffe(block_size*(i-1)+1:block_size*i));
    
    %% lasso polynomial evaluation
    D_t=x2fx(tBlock_var,'linear');
    coef=coef_saver(:,i); 
    coef0=coef0_saver(:,i);
    Z_appro_part1=D_t*coef+coef0;
    

    
%     indrow=indrow_saver(:,i);
    BigEval=tBlock_var';
    % sparsity selection
    [indrow,~]=find(coef>0);
    BigEval=BigEval(indrow-1,:);
    weights_tmp=weights_saver(:,i);
    
%     H_len=H_len+indrow;
    H_tmp=H_saver(H_len+1:H_len+length(indrow),:);
    H_len=H_len+length(indrow);
    Z_appro_part2=weights_tmp'*exp(-1/r*(repmat(diag(H_tmp'*H_tmp),1,size(BigEval,2))-2*H_tmp'*BigEval+repmat(diag(BigEval'*BigEval)',batch_size,1)));
    
    Z_appro_part=Z_appro_part1+Z_appro_part2';
    Z_appro=Z_appro+Z_appro_part;
    Z_appro_part_saver=[Z_appro_part_saver,Z_appro_part];
end

%% Relative error
Re=sum(abs(Z_appro-M_Y_test))/length(M_Y_test);
%% full connect network
% reforce

% bias=mean(M_Y_test-Z_appro);
% Z_appro=Z_appro+bias;

% LS
% Z_appro=coef_LS(1)+coef_LS(2)*Z_appro_part_saver(:,1)+coef_LS(3)*Z_appro_part_saver(:,2);

% MLP
testPredictions = predict(Mdl,Z_appro_part_saver);
% plot(YTest,testPredictions,".")
% hold on
% plot(YTest,YTest)
% hold off
% xlabel("True MPG")
% ylabel("Predicted MPG")
%% Relative error calc
RE1=sum((abs(M_Y_test-testPredictions)))/length(M_Y_test);
RE2=immse(M_Y_test,testPredictions);
