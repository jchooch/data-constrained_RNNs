%%% train_100_sparse_nets.m

close all;
clear all;
clc;

load('calact410Nstd.mat')
load('inputs.mat')

% find optimal lambda via hyperparameter search

lambda = 0;

for net_index = 1:100
    [J, R, N, T, varData, chi2, pVars] = dcRNN_train(calact410Nstd, inputs, 50, lambda);
    network.J = J; network.R = R; network.N = N; network.T = T;
    network.varData = varData; network.chi2 = chi2; network.pVars = pVars;
    save(['sparse_ensemble_20220713/sparse_output_', num2str(net_index), '.mat'], 'network');
    clear network
end
