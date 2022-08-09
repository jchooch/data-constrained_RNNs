%%% analyse_error_ensemble.m

clear all;
close all;
clc;

load("mode_5_rep_1.mat");
J_1 = network.J;
load("mode_5_rep_2.mat");
J_2 = network.J;
load("mode_5_rep_3.mat");
J_3 = network.J;

load("mode_1_rep_1.mat");
J_4 = network.J;
load("mode_1_rep_2.mat");
J_5 = network.J;
load("mode_1_rep_3.mat");
J_6 = network.J;

test_Js = {J_1, J_2, J_3};
vanilla_Js = {J_4, J_5, J_6};

figure(1);
for rep = 1:3
    subplot(1, 3, rep);
    h = heatmap(test_Js{rep});
    caxis([-1, 1]);
    h.Colormap = redblue(100);
    h.GridVisible = 'off';
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    h.Title = sprintf('Mode 4 Rep %i', rep);
    h.XLabel = 'Presynaptic';
    h.YLabel = 'Postsynaptic';
end

figure(2);
for rep = 1:3
    subplot(1, 3, rep);
    h = heatmap(vanilla_Js{rep});
    caxis([-1, 1]);
    h.Colormap = redblue(100);
    h.GridVisible = 'off';
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    h.Title = sprintf('Mode 1 Rep %i', rep);
    h.XLabel = 'Presynaptic';
    h.YLabel = 'Postsynaptic';
end

figure(3);
test_weights = [J_1; J_2; J_3];
vanilla_weights = [J_4; J_5; J_6];
title(sprintf('Weight Distribution for Different L1 Lambdas'));
hold on;
histogram(vanilla_weights, 'BinEdges', linspace(-1, 1, 200), 'FaceAlpha', 0.3, 'FaceColor', 'b')
histogram(test_weights, 'BinEdges', linspace(-1, 1, 200), 'FaceAlpha', 0.3, 'FaceColor', 'r')
%ylim([0 5000]);
legend('No L1 Penalty', 'Sensitive L1 Penalty')
hold off;