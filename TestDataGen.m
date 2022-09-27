%% intial
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

%% LHS Sampling
Num=500; Dim=6;

TMP_RAND=-4+8*lhsdesign(Num,Dim);
M_X_test=repmat(mu,Num,1)+repmat(sigma,Num,1).*TMP_RAND;
M_Y_test=run_file(analysis.hspicepath,M_X_test,Num);
save('M_test','M_X_test','M_Y_test');
%% MC Sampling
Num=500; Dim=6;

TMP_RAND=randn(Num,Dim);
M_X_test_MC=repmat(mu,Num,1)+repmat(sigma,Num,1).*TMP_RAND;
M_Y_test_MC=run_file(analysis.hspicepath,M_X_test_MC,Num);
save('M_test_MC','M_X_test_MC','M_Y_test_MC');

%% Truncated MC sampling
Num=500; Dim=6;
pretruncMean=0; pretruncSD=1;
untruncated = makedist('Normal',pretruncMean,pretruncSD);
truncated = truncate(untruncated,pretruncMean-3,pretruncMean+3);
TMP_CRAND=random(truncated,Num,Dim);
M_X_test_TMC=repmat(mu,Num,1)+repmat(sigma,Num,1).*TMP_CRAND;
M_Y_test_TMC=run_file(analysis.hspicepath,M_X_test_TMC,Num);
save('M_test_TMC','M_X_test_TMC','M_Y_test_TMC');


%############################### For transfer learning ###################%
%% Truncated MC sampling with different voltage
%% 0p8
Num=500; Dim=6;
pretruncMean=0; pretruncSD=1;
untruncated = makedist('Normal',pretruncMean,pretruncSD);
truncated = truncate(untruncated,pretruncMean-3,pretruncMean+3);
TMP_CRAND=random(truncated,Num,Dim);
M_X_test_TMC0p8=repmat(mu,Num,1)+repmat(sigma,Num,1).*TMP_CRAND;
M_Y_test_TMC0p8=run_file(analysis.hspicepath,M_X_test_TMC0p8,Num);
save('M_test_TMC0p8','M_X_test_TMC0p8','M_Y_test_TMC0p8');

%% 0p6
% please change the voltage of netlist artificially
M_X_test_TMC0p6=M_X_test_TMC0p8;
M_Y_test_TMC0p6=run_file(analysis.hspicepath,M_X_test_TMC0p8,Num);
save('M_test_TMC0p6','M_X_test_TMC0p6','M_Y_test_TMC0p6');
%% 1p0
% please change the voltage of netlist artificially
M_X_test_TMC1p0=M_X_test_TMC0p8;
M_Y_test_TMC1p0=run_file(analysis.hspicepath,M_X_test_TMC0p8,Num);
save('M_test_TMC1p0','M_X_test_TMC1p0','M_Y_test_TMC1p0');

