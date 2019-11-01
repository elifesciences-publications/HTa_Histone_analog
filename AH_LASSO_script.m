%Loading and formating the input table (an example is provided): 

load 'Normalized_MNase_and_predictors_Tacidophilum_40_65bp_TEAC_r2d2_80_predictors.mat'

inputTable=struct2table(NormInputTrainExport);
predictorNames_prep = fieldnames(NormInputTrainExport);

predictorNames=predictorNames_prep(1:length(predictorNames_prep)-1);
predictors = table2array(inputTable(:, predictorNames));
response = inputTable.MNaseData;

%Applying the built in lasso algorithm, 10 xinternal cross validation
[B,FitInfo] = lasso(predictors,response,'CV',10,'PredictorNames',predictorNames);

%Plot the cross-validated fits.
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend
savefig('MSE_CV_Normalized_MNase_and_predictors_Tacidophilum_40_65bp_TEAC_r2d2_80_predictors.fig')
idxLambda1SE = FitInfo.Index1SE;
coef = B(:,idxLambda1SE);
coef0 = FitInfo.Intercept(idxLambda1SE);

%Export the coefficients for later use
T=table(predictorNames,coef)
T.Properties.VariableNames={'Name','Value'}
T2=table('Intercept',coef0)
T2.Properties.VariableNames={'Name','Value'}

T3=vertcat(T,T2)

%Output : Coefficient that will be used to build predicted track 
writetable(T3,'Coefficients_Normalized_MNase_and_predictors_Tacidophilum_40_65bp_TEAC_r2d2_80_predictors.txt','Delimiter','\t') 

%Output : All fitting informations
writetable(struct2table(FitInfo),'FitInfo_Normalized_MNase_and_predictors_Tacidophilum_40_65bp_TEAC_r2d2_80_predictors.txt','Delimiter','\t') 