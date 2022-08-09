%%% train_lambda_ensemble_ver1.m

clc;
close all;
clear all;

load('calact410Nstd.mat');
load('inputs.mat');
load('lambdas.mat');

repetitions = 10; % number of networks trained per lambda
Js = {}; Rs = {}; chi2s = {}; pVarss = {};
chi2_means = []; pVars_means = []; % statistics over reptitions per lambda

chi2_sum = 0;
pVars_sum = 0;

index_row = 1; % lambda
for lambda = lambdas
    fprintf('Lambda = %.10g\n', lambda)
    index_col = 1; % repetition
    for repetition = 1:repetitions
        [J, R, N, T, varData, chi2, pVars] = dcRNN_train_L1_ver1(calact410Nstd, inputs, 100, lambda);
        Js{index_row, index_col} = J;
        Rs{index_row, index_col} = R;
        chi2s{index_row, index_col} = chi2;
        pVarss{index_row, index_col} = pVars;
        index_col = index_col + 1;
        chi2_sum = chi2_sum + chi2(end);
        pVars_sum = pVars_sum + pVars(end);
    end
    index_row = index_row + 1;
    chi2_mean = chi2_sum / repetitions;
    pVars_mean = pVars_sum / repetitions;
    chi2_means = [chi2_means; chi2_mean];
    pVars_means = [pVars_means; pVars_mean];
end

save('20_lambs_10_reps_L1_ver1.mat')