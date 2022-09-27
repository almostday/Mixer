clc;clear;

analysis.path='D:\MATLAB\R2021a\workspace\BCD';
analysis.hspicepath = '"C:\synopsys\Hspice_B-2007.03\BIN\hspice.exe"';
cd(analysis.path);
% parameters set
[analysis.model.mean,analysis.model.sigma]=initial_format();
mu=analysis.model.mean;
sigma=analysis.model.sigma;

% model_save={ };
% 选取切割点并计算其响应
    % 初始点选用所有维度中心位置的近似值
len=length(mu);

cut_point_norm=zeros(1,6);
cut_point=cut_point_norm.*sigma+mu;
ft=run_file(analysis.hspicepath,cut_point,1);

X_train=[];
Y_train=[];
%% single variable set
% para_save=zeros(len,10);
% LHS sampling
S_X_train=[]; S_Y_train=[];
for i=1:len  
     lhs_num=300;
     nor_num=200;
     num=lhs_num+nor_num;
     tmp3_norm=-4+4*lhsdesign(num,1);
     tmp3=tmp3_norm.*repmat(sigma(i),num,1)+repmat(mu(i),num,1);
     sample=repmat(cut_point,num,1);
     sample(:,i)=tmp3;
     y_tmp3=run_file(analysis.hspicepath,sample,num);
     S_X_train=[S_X_train;tmp3_norm];
     S_Y_train=[S_Y_train;y_tmp3];
end
S_X_train=reshape(S_X_train,[num,len]);
S_Y_train=reshape(S_Y_train,[num,len]);
%% multi-variable set
Num=1000; Dim=6;

TMP_RAND=-6+12*lhsdesign(Num,Dim);
M_X_train=repmat(mu,Num,1)+repmat(sigma,Num,1).*TMP_RAND;
M_Y_train=run_file(analysis.hspicepath,M_X_train,Num);

%% failure center
% load('ymy01_failed_sample_norm','failed_sample_norm');
% load('failed_sample_norm0p9','failed_sample_norm0p9');
% failed_sample_norm=failed_sample_norm0p9;
% cut_mean=mean(failed_sample_norm);
% cut_point_norm=cut_mean;
% cut_point_norm=failed_sample_norm(1,:);
% % 
% lower=min(failed_sample_norm);
% upper=max(failed_sample_norm);
% [t1,t2]=size(failed_sample_norm);
% for t=1:t1
%      failed_sample(t,:)=mu+failed_sample_norm(t,:).*sigma;
% end

