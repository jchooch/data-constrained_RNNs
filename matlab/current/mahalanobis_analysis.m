%%% mahalanobis_analysis.m

clear all;
close all;
clc;

load('stim_on_or_off.mat');
load('reduced_currents.mat');
load('inputs.mat');

inputs = logical(inputs);

[binary_starts, binary_stops] = groupLims(groupTrue(stim_on_or_off), 1);
starts_and_stops = cell(1, size(inputs, 1));

for channel = 1:size(inputs, 1)
    [starts, stops] = groupLims(groupTrue(inputs(channel, :)), 1);
    starts_and_stops{1,channel} = [starts; stops];
end

mahal_distances = cell(size(reduced_currents));
mean_currents = cell(size(reduced_currents));

for iT = 1:size(reduced_currents, 1)
    for iS = 1:size(reduced_currents, 2)
        mean_currents{iT,iS} = mean(reduced_currents{iT,iS}, 2);
        for t = 1:length(reduced_currents{iT,iS})
            mahal_distances{iT,iS}(t) = mahalanobis_distance(reduced_currents{iT,iS}(:,t), mean_currents{iT,iS});
        end
    end
end

patch_colors = [0.6549, 0.7804, 0.9059;
                    1, 0.8941, 1;
                    1, 0.93, 0.85;
                    0.8745, 0.9824, 0.8863;
                    1, 0.85, 0.82;
                    0.72, 0.92, 0.97;
                    0.9020, 0.9020, 0.9804;
                    1, 0.9725, 0.9059];

channel_names = ["REUM", "RERM", "REDM", "RELM", "LEUM", "LERM", "LEDM", "LELM"];

%{
figure()
hold on;
patch([0 1 1 0], [1 1 0 0], patch_colors(1,:))
patch([1 2 2 1], [1 1 0 0], patch_colors(2,:))
patch([2 3 3 2], [1 1 0 0], patch_colors(3,:))
patch([3 4 4 3], [1 1 0 0], patch_colors(4,:))
patch([4 5 5 4], [1 1 0 0], patch_colors(5,:))
patch([5 6 6 5], [1 1 0 0], patch_colors(6,:))
patch([6 7 7 6], [1 1 0 0], patch_colors(7,:))
patch([7 8 8 7], [1 1 0 0], patch_colors(8,:))
hold off;
%}

patch_mode = 2; % 1=binary(on/off), 2=all
figure(1)
hold on
%{
for iT = 1:size(reduced_currents, 1)
    for iS = 1:size(reduced_currents, 2)
        plot(linspace(1, length(mahal_distances{iT,iS}), length(mahal_distances{iT,iS})), mahal_distances{iT,iS})
    end
end
%}
newcolors = [0, 0.4470, 0.7410; 0.8549, 0.1725, 0.2627;
                0.3010, 0.7450, 0.9330; 1.0000, 0.3882, 0.2784];
colororder(newcolors);
for iT = 1:size(reduced_currents, 1)
    for iS = 1:size(reduced_currents, 2)
        plot(linspace(1, length(mahal_distances{iT,iS}), length(mahal_distances{iT,iS})), mahal_distances{iT,iS}, 'LineWidth', 2)
    end
end
if patch_mode == 1
    for i = 1:length(binary_starts)
        patch([binary_starts(i) binary_stops(i) binary_stops(i) binary_starts(i)], [max(ylim) max(ylim) 0 0], [0.9290 0.6940 0.1250])
    end
elseif patch_mode == 2
    for channel = 1:size(inputs, 1)
        these_starts = starts_and_stops{1,channel}(1,:);
        these_stops = starts_and_stops{1,channel}(2,:);
        for i = 1:length(these_starts)
            patch([these_starts(i) these_stops(i) these_stops(i) these_starts(i)], [max(ylim) max(ylim) 0 0], patch_colors(channel,:))
            text(these_starts(i)+3, 30, channel_names(channel), 'Color', 'k')
        end
    end
end
for iT = 1:size(reduced_currents, 1)
    for iS = 1:size(reduced_currents, 2)
        plot(linspace(1, length(mahal_distances{iT,iS}), length(mahal_distances{iT,iS})), mahal_distances{iT,iS}, 'LineWidth', 2)
    end
end
xlabel('Timepoint')
ylabel('Mahalanobis distance')
legend('Leftward to Leftward', 'Rightward to Leftward', 'Leftward to Rightward', 'Rightward to Rightward')
%legend('Leftward to Leftward', 'Rightward to Leftward', 'Leftward to Rightward', 'Rightward to Rightward', ...
    %'REUM', '', '', 'RERM', '', '', '', '', '', '', '', '', '', 'REDM', '', '', '', 'RELM', 'LEUM', 'LERM', 'LEDM', 'LELM')
xlim([1 1000])
hold off

% "REUM", "RERM", "REDM", "RELM", "LEUM", "LERM", "LEDM", "LELM"





