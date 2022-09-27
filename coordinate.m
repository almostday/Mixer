rng default
X = randn(100,10);
weights = [0.6;0.5;0.7;0.4];
y = X(:,[2 4 5 7])*weights + randn(100,1)*0.1; % Small added noise

B = lassoglm(X,y);