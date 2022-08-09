% see_activity_fit.m

function see_activity_fit(R, Z, neurons)
    disp('Size of RNN activities:')
    disp(size(R))
    disp('Size of ZNN activities:')
    disp(size(Z))
    T = size(R, 2);
    disp('Neurons:')
    disp(neurons)
    figure('Position',[100 100 1800 600])
    for neuron_index = 1:length(neurons)
        subplot(length(neurons), 1, neuron_index)
        plot(linspace(1, T, size(Z, 2)), Z(neuron_index,:), '-b', 'DisplayName', 'znn');
        hold on;
        plot(linspace(1, T, size(R, 2)), R(neuron_index,:), '-r', 'DisplayName', 'rnn');
        hold off;
        lgd = legend;
        lgd.NumColumns = 1;
        lgd.Location = 'northeast';
        title(sprintf('ZNN vs RNN Neuron %i', neurons(neuron_index)));
        xlabel('Time');
        ylabel('Activity / \DeltaF/F');
        set(gca,'Box','off','TickDir','out','FontSize',9);
        %{
        subplot(3,6,[13 14 15 16 17 18], 'align');
        h = heatmap(tStim(:,1:2:end));
        grid off
        axs = struct(gca); %ignore warning that this should be avoided
        cb = axs.Colorbar;
        cb.Ticks = 0:4;
        cb.TickLabels = {'No stim','Right','Forward','Left','Backward'};
        XLabels = 1:1000;
        CustomXLabels = string(XLabels);
        CustomXLabels(mod(XLabels,100) ~= 0) = " ";
        h.XDisplayLabels = CustomXLabels;
        h.YData = {'Left Eye','Right Eye'};
        mymap = [1 1 1
            1 0 0
            0 1 0
            0 0 1
            1 0 1];
        colormap(mymap);
        colorbar;
        axs.Fontsize = 10;
        drawnow;
        saveas(f,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
        %}
    end
end

%{
f7 = figure('Position',[100 100 1800 600]);
clf(f7);
figure(7);
Jhm = heatmap(J);
S = struct(Jhm); % Undocumented
ax = S.Axes;    % Undocumented
% Remove grids
Jhm.GridVisible = 'off';
% Place lines around selected columns and row
% Assumes columns and rows are 1 unit in size!
col = [54 205 271];    
row = [54 205 271];
arrayfun(@(x)xline(ax,x,'k-','Alpha',1),[col+.5]);
arrayfun(@(x)yline(ax,x,'k-','Alpha',1),[row+.5]);
%xline(ax, [col-.5,col+.5], 'k-', 'Alpha'); % see footnotes [1,2]
%yline(ax, [row-.5,row+.5], 'k-', 'Alpha'); % see footnotes [1,2]
%grid on
axsJ = struct(gca); %ignore warning that this should be avoided
%cb = axs.Colorbar;
%cb.Ticks = 0:4;
%cb.TickLabels = {'No stim','Right','Forward','Left','Backward'};
XLabels = 1:410;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,500) ~= 0) = " ";
CustomXLabels(XLabels == 27) = "Forward";
CustomXLabels(XLabels == 130) = "Right";
CustomXLabels(XLabels == 239) = "Backward";
CustomXLabels(XLabels == 344) = "Left";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
Jhm.XDisplayLabels = CustomXLabels;
YLabels = 1:410;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,500) ~= 0) = " ";
CustomYLabels(YLabels == 27) = "Forward";
CustomYLabels(YLabels == 130) = "Right";
CustomYLabels(YLabels == 239) = "Backward";
CustomYLabels(YLabels == 344) = "Left";
Jhm.YDisplayLabels = CustomYLabels;
colormap(cmap);
colorbar;
caxis([-.65 .65]);
axs.Fontsize = 10;
Jhm.Title = 'Synaptic Strength Weights';
drawnow;
saveas(f7,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
%}