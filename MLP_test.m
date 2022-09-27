load carbig
X = [Acceleration Cylinders Displacement Weight];
Y = MPG;

%%
R = rmmissing([X Y]);
X = R(:,1:end-1);
Y = R(:,end);

%%
rng("default") % For reproducibility of the partition
c = cvpartition(length(Y),"Holdout",0.20);
trainingIdx = training(c); % Indices for the training set
XTrain = X(trainingIdx,:);
YTrain = Y(trainingIdx);
testIdx = test(c); % Indices for the test set
XTest = X(testIdx,:);
YTest = Y(testIdx);
%%

Mdl = fitrnet(XTrain,YTrain,"Standardize",true, ...
    "LayerSizes",[30 10]);

%%
testMSE = loss(Mdl,XTest,YTest);
testPredictions = predict(Mdl,XTest);
plot(YTest,testPredictions,".")
hold on
plot(YTest,YTest)
hold off
xlabel("True MPG")
ylabel("Predicted MPG")


