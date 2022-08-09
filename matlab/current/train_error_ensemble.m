%%% train_error_ensemble.m
%%% RENAME TO train_L1_ensemble.m FOR CLARITY

clear all;
close all;
clc;

load('inputs.mat');
load('calact410Nstd.mat');
lambdas = [1e-20, 1e-15, 1e-10, 1e-5];

for error_mode = [2, 3, 4, 6, 7]
    for rep = 1:50
        disp('Error mode, repetition: ')
        disp([error_mode, rep])
        [J, R, N, T, varData, chi2, pVars] = dcRNN_train(calact410Nstd, inputs, 50, error_mode, 1e-10);
        network.J = J; network.R = R; network.N = N; network.T = T;
        network.varData = varData; network.chi2 = chi2; network.pVars = pVars;
        save(sprintf('error_ensemble/mode_%i_rep_%i.mat', error_mode, rep), 'network');
        clear network;
    end
end

for error_mode = [1, 5]
    for lambda = lambdas
        for rep = 1:50
            disp('Error mode, repetition, lambda: ')
            disp([error_mode, rep, lambda])
            [J, R, N, T, varData, chi2, pVars] = dcRNN_train(calact410Nstd, inputs, 50, error_mode, lambda);
            network.J = J; network.R = R; network.N = N; network.T = T;
            network.varData = varData; network.chi2 = chi2; network.pVars = pVars;
            save(sprintf('error_ensemble/mode_%i_lamb_%.0e_rep_%i.mat', error_mode, lambda, rep), 'network');
            clear network;
        end
    end
end

