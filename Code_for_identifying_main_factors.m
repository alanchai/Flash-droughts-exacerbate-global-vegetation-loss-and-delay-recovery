clc;
clear; 
% Load data files
%load datax.mat;
%load datay.mat;


t = templateTree('NumVariablesToSample', 'all', ...
    'PredictorSelection', 'interaction-curvature', 'Surrogate', 'on');
Mdl3 = fitrensemble(datax, datay, 'Method', 'Bag', 'NumLearningCycles', 500, ...
    'Learners', t);
impOOB = oobPermutedPredictorImportance(Mdl3);
figure(2)
bar(impOOB);
title('Unbiased Predictor Importance Estimates');
xlabel('Predictor variable');
ylabel('Importance');
h = gca;
h.XTickLabel = Mdl3.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
trees = 500;
leaf = 5;
OOBPrediction = 'on';
OOBPredictorImportance = 'on';
Method = 'regression';
Mdl4 = TreeBagger(trees, datax, datay, ...
    'OOBPredictorImportance', 'on', ...
    'Method', Method, ...
    'OOBPrediction', OOBPrediction, ...
    'MinLeaf', leaf);
importance = Mdl4.OOBPermutedPredictorDeltaError;
figure(3)
bar(importance);
title('OOB Permuted Predictor Importance Estimates'); 
figure(1)
for i = 1:9
    subplot(3, 3, i);
    plotPartialDependence(Mdl4, i); 
    title(['Partial Dependence for Predictor ' num2str(i)]);
end