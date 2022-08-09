%%% compare_normal_and_sparse_Js.m

clear all;
close all;
clc;

load('ensemble/output_1.mat');
J_1 = J;
load('ensemble/output_50.mat');
J_2 = J;
load('20_lambs_10_reps_L1_ver1.mat')

figure(1); 
subplot(2,2,1);
h1 = heatmap(Js{8, 1}); 
h1.Colormap = redblue(100); 
h1.GridVisible = 'off'; 
caxis([-1 1]); 
Ax = gca; 
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h1.Title = 'Sparse J (Lambda = 1e-13)';
h1.XLabel = 'Presynaptic';
h1.YLabel = 'Postsynaptic';

subplot(2,2,2); 
h2 = heatmap(J_1); 
h2.Colormap = redblue(100); 
h2.GridVisible = 'off'; 
caxis([-1 1]); 
Ax = gca; 
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h2.Title = 'Normal J (J_1)';
h2.XLabel = 'Presynaptic';
h2.YLabel = 'Postsynaptic';

subplot(2,2,3);
h3 = heatmap(Js{12, 1}); 
h3.Colormap = redblue(100); 
h3.GridVisible = 'off'; 
caxis([-1 1]); 
Ax = gca; 
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h3.Title = 'Sparse J (Lambda = 1e-9)';
h3.XLabel = 'Presynaptic';
h3.YLabel = 'Postsynaptic';

subplot(2,2,4); 
h4 = heatmap(J_2); 
h4.Colormap = redblue(100); 
h4.GridVisible = 'off'; 
caxis([-1 1]); 
Ax = gca; 
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h4.Title = 'Normal J (J_2)';
h4.XLabel = 'Presynaptic';
h4.YLabel = 'Postsynaptic';

figure(2)
title(sprintf('Weight Distribution for Different L1 Lambdas'));
hold on;
%histogram(Js{1,1}, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.2, 'FaceColor', 'r')
%histogram(Js{5,1}, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.3, 'FaceColor', 'g')
histogram(Js{8,1}, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.2, 'FaceColor', 'b')
%histogram(J_1, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.2, 'FaceColor', 'yellow')
histogram(J_2, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.2, 'FaceColor', 'g')
%histogram(Js{15,1}, 'BinEdges', linspace(-1, 1, 500), 'FaceAlpha', 0.3, 'FaceColor', 'y')
ylim([0 5000]);
legend('Lambda = 1e-13', 'Normal J 2')
hold off;

figure(3)
for lamb = 1:20
    for rep = 1:10
        ...
    end
end
means = [mean(abs(Js{8, 1}), 'all'), mean(abs(Js{12, 1}), 'all'), mean(abs(J_1), 'all'), mean(abs(J_2), 'all')];

bar(means); 

set(gca, 'XTickLabel', ["1e-13" "1e-9" "0" "0"]); xlabel('Lambda');



