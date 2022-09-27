load acetylene


X = [x1 x2 x3];
D = x2fx(X,'interaction');
% D(:,1) = []; % No constant term
rng default % For reproducibility 
[B,FitInfo] = lasso(D,y,'CV',10);

