clc;
clear all;
close all;

RDM_srce_past=randi([1 1000],[100 1])*0.001;
RDM_targ_past=randi([1 1000],[100 1])*0.001;
RDM_targ_prsnt=randi([1 1000],[100 1])*0.001;

%% Linear regression
[~, goodness_without, ~] = fit(RDM_targ_past,RDM_targ_prsnt,'poly1');
[~, goodness_with, ~] = fit([RDM_targ_past,RDM_srce_past],RDM_targ_prsnt,'poly12');
unexplained_variance=log((1-goodness_without.adjrsquare)./(1-goodness_with.adjrsquare))

%% General Linear Model

T = table(RDM_targ_prsnt,RDM_targ_past,RDM_srce_past);
T.Properties.VariableNames = {'RDM_targ_prsnt' 'RDM_targ_past' 'RDM_srce_past'};
      

glme_without = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_srce_past','Link','identity','DummyVarCoding','effects');
glme_with = fitglme(T,'RDM_targ_prsnt ~ 1 + RDM_targ_past + RDM_srce_past','Link','identity','DummyVarCoding','effects')

unexplained_variance=log((1-glme_without.Rsquared.Adjusted)./(1-glme_with.Rsquared.Adjusted))
