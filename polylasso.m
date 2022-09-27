function [coef,coef0] = polylasso(traindata,Z_samples)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明


[B,FitInfo] = lasso(traindata,Z_samples,'CV',10);
idxLambda1SE = FitInfo.Index1SE;
coef = B(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);
end

