%% Load validation dataset
load CPall_CPhrs.mat

%% Observed phase
y=dat.time;

%% Load the model
load model_NC5_G100_1sample_CPhrs.mat

%% z-score within samples
dat.num=zscore(dat.num);
[n,p]=size(dat.num);

%% Predict using the models parameters for the linear regression
yfit = [ones(n,1) dat.num(:,model.index(1:model.max_var))]*model.beta;
yfit = postprocess_sincos(yfit);
rsquared = evaluate_perf(y,yfit)

%% Save results
xlswrite('Prediction_PLSR_1sample_5_100_CPhrs.csv',[y, yfit],'Observed Predicted');


