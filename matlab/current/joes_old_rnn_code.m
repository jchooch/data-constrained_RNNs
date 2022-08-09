function  joes_old_rnn_code()
    plotStatus = 1;
    simfilename = '/Users/joechoo-choy/Documents/github_jchooch/zfish_pretectum_rnn/matlab/current/calact410Nstd.mat'; %replace with input path
    outputname = sprintf('/Users/joechoo-choy/Documents/github_jchooch/zfish_pretectum_rnn/output_last');
    %outputname ='/Users/joechoo-choy/Downloads/output_12';
    fname = '/Users/joechoo-choy/Downloads/';
    %colmap = '/Users/joechoo-choy/Downloads/redblue/redblue.mat';
    learn = 1; %lets RNN learn
    nRunTot = 50; %number of steps
    nFree = 5; %how many epochs to freeze learning before final epoch
    nFreePre = 0; %how many epochs from first one that learning is frozen
    firing_rates = load(simfilename); %hold-over name - call it whatever - it gets converted into Adata later
    firing_rates = firing_rates.calact410Nstd;
    %cmap = load('/Users/joechoo-choy/Downloads/redblue/redblue.mat');
    %cmap = redblue(100);
    %cmap = darkb2r(-3,3);
    if (exist('nFreePre') == 0)
      nFreePre = 0;
    end
    %{
    function c = redblue(m)
        if nargin < 1, m = size(get(gcf,'colormap'),1); end
        if (mod(m,2) == 0)
            % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = m*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
        else
            % From [0 0 1] to [1 1 1] to [1 0 0];
            m1 = floor(m*0.5);
            r = (0:m1-1)'/max(m1,1);
            g = r;
            r = [r; ones(m1+1,1)];
            g = [g; 1; flipud(g)];
            b = flipud(r);
        end
        c = [r g b];
    end
    %}

    cmap = redblue(100);

    %variable values
    
    % control parameter - if positive, ntwk produces spontaneous activity
    % with non-trivial dynamics. if negative, it does not (CS)
    g = 1.25;   % Andalman used between 1.2 and 1.3 for g. Generally recommended to be between 1 and 1.5
    dtData = 0.5; % ? frame rate given in Supplemental Info (CS) - 1/frame rate - our actual frame rate is slightly faster than 2 so actual dtData is somewhat lower
    dt = 0.25; % integration step in sec - this value was used in Andalman
    tau = 1; % 2.5s time constant of each unit in the network - value was used in Andalman - analogous to neuron decay time
    P0 = 1; % free parameter to the learning rule for J (CS) - Recommended to be 1 to 0.01 from Sussillo & Abbott
    N = 410; % #neurons - however many rows or neuron time series you are training
    tauWN = 1; % correlation time of external inputs to network in sec (CS)
    ampIn = 0.005; % scale of external inputs to network (h_0) (CS)

    selectivities = zeros(1, N);
    selectivities(1:54) = 1; % forward binocular
    selectivities(55:118) = 2; % right lateral
    selectivities(119:205) = 3; % right medial
    selectivities(206:268) = 4; % backward binocular
    selectivities(269:350) = 5; % left medial
    selectivities(351:410) = 6; % left lateral
    
    %Frames and Time length
    %if dtData = 0.5 then data interval is twice length of start:end.
    data_start_0ndx = 0;
    data_end_0ndx = 1999;
    data_start_time = 0;
    data_end_time = 1000;
    
    N = double(N);

    tData_ind = (data_start_0ndx+1):(data_end_0ndx+1);  
    tData = [0:dtData:1999*dtData]; %How many times J is updated
    
    t = data_start_time:dt:data_end_time;  %Number of integration time steps
    
    lastchi2s = zeros(length(tData), 1);

    iModelSample = zeros(length(tData), 1);
    for i=1:length(tData)
        [tModelSample, iModelSample(i)] = min(abs(tData(i)-t));
    end    
    
    Adata = firing_rates; %Real calcium data used for initial conditions and for updating J in err.
    
    Adata = Adata/max(max(Adata));
    
    Adata = min(Adata, 0.999);
    Adata = max(Adata, -0.999);
    
    stdData = std(reshape(Adata, N*length(tData_ind), 1));
    
    %Gaussian white noise with amplitude no more than 1% of baseline
    %according to Andalman
    %set seed - have random number gen and then randn
    ampWN = sqrt(tauWN/dt);
    iWN = ampWN*randn(N, length(t));
    inputN = ones(N, length(t));
    
    for tt = 2: length(t)
        inputN(:, tt) = iWN(:, tt) + (inputN(:, tt - 1) - iWN(:, tt))*exp(-(dt/tauWN)); %white noise
    end
    inputN = ampIn*inputN;
    
    %Initialize 
    J0 = g*randn(N,N)/sqrt(N); %synaptic strength matrix
    J = J0;
    R = NaN(N, length(t)); %nonlinear activity
    AR = NaN(N+8, length(t)); %augmented activity for training input wts
    JR = zeros(N, 1); %product of synaptic strength weights and nonlinear activity
    %inputs - e.g. LU means Left eye Up stim
    [inputLR,inputLU,inputLL,inputLD,inputRR,inputRU,inputRL,inputRD] = deal(zeros(N,length(t)));
    
    %external input weights - make sure using same seed for noise and
    %external input
    
    %CHANGE TO JUST wLR not wLR0
    [wLR0,wLU0,wLL0,wLD0,wRR0,wRU0,wRL0,wRD0] = deal(randn(N,1)/sqrt(N));
    wLR = wLR0;wLU = wLU0;wLL = wLL0;wLD = wLD0;wRR = wRR0;wRU = wRU0;wRL = wRL0;wRD = wRD0;
    
    %stim on - manually coded times for now since there was some slight
    %difference in times between and during stimuli in recording.
    %Generally, stimuli last ~16 frames or 8 real seconds or 32 timesteps
    [tLR,tLU,tLL,tLD,tRR,tRU,tRL,tRD] = deal(zeros(1,length(t)));
    tLR(1184:1214)=1; tLR(1752:1782)=1; tLR(2774:2804)=1;
    tLU(162:190)=1; tLU(502:532)=1; tLU(730:760)=1; tLU(1070:1100)=1; tLU(1638:1668)=1; tLU(2206:2236)=1; tLU(2548:2578)=1; tLU(2662:2690)=1; tLU(2888:2918)=1; tLU(3116:3146)=1;
    tLL(48:78)=1;tLL(1980:2010)=1;tLL(2320:2350)=1;tLL(3002:3032)=1;
    tLD(274:304)=1;tLD(842:872)=1;tLD(1412:1440)=1;tLD(2434:2464)=1;tLD(3230:3260)=1;tLD(3684:3714)=1;
    tRR(1184:1214)=1;tRR(1752:1782)=1;tRR(2774:2804)=1;
    tRU(502:532)=1;tRU(616:646)=1;tRU(956:986)=1;tRU(1070:1100)=1;tRU(1298:1328)=1;tRU(1412:1440)=1;tRU(1524:1554)=1;tRU(2548:2578)=1;tRU(2662:2690)=1;
    tRL(48:78)=1;tRL(1980:2010)=1;tRL(2320:2350)=1;tRL(3002:3032)=1;
    tRD(274:304)=1;tRD(730:760)=1;tRD(1638:1668)=1;tRD(1866:1896)=1;tRD(2092:2122)=1;tRD(2206:2236)=1;tRD(3116:3146)=1;tRD(3342:3372)=1;tRD(3570:3600)=1;tRD(3798:3828)=1;tRD(3912:3940)=1;
    
    t_hm = [0:0.5:1999*0.5];
    load('left_stim_times.mat')
    load('right_stim_times.mat')
    tStim = [left_stim_times; right_stim_times];
    
    %amp - Appears Andalman et al. recommends 10x larger than amplitude of
    %noise
    amp_rgc = 0.05;
    
    if (learn)
        PJ = P0 * eye(N+8, N+8);
    end
    
    [data_units, data_steps] = size(Adata);
    fprintf('data_steps: %d\n', data_steps);
    fprintf('dt: %f\n', dt);
    fprintf('length(t): %d\n', length(t));
    assert(data_steps >= length(t) * dt);

    chi2 = zeros(1, nRunTot);
    pVars = zeros(1, nRunTot);

    for nRun = 1:nRunTot %Number of epochs
        H = Adata(:, 1); %Initialize activities of all neurons with real values at first time point
        R(:, 1) = tanh(H); %nonlinearly transformed activities
        tLearn = 0; %param for when to update J matrix
        iLearn = 1; %Used to index Adata to subtract from model predicted rates in err.
        [epoch_LR,epoch_LU,epoch_LL,epoch_LD,epoch_RR,epoch_RU,epoch_RL,epoch_RD] = deal(0); %set epochs of external inputs to 0
     
        for tt = 2:length(t)-2 %time steps for each epoch. 
            tLearn = tLearn + dt; %update for each time step
            R(:, tt) = tanh(H); %nonlinear transformation of activities
            if tLR(tt) == 1 %conditionals for training external input weights when stim is on
              epoch_LR = 1;
            end
            if tLU(tt) == 1
              epoch_LU = 1;
            end
            if tLL(tt) == 1
              epoch_LL = 1;
            end
            if tLD(tt) == 1
              epoch_LD = 1;
            end
            if tRR(tt) == 1
              epoch_RR = 1;
            end
            if tRU(tt) == 1
              epoch_RU = 1;
            end
            if tRL(tt) == 1
              epoch_RL = 1;
            end
            if tRD(tt) == 1
              epoch_RD = 1;
            end
            
            JR = J*R(:, tt) + inputN(:, tt); %product of synaptic strength weight matrix and nonlinear activities + gaussian noise
            
            %external inputs at each time step - input on if time vectors
            %are 1 at time tt.
            inputLR(:, tt) = amp_rgc*tLR(tt).*wLR;
            inputLU(:, tt) = amp_rgc*tLU(tt).*wLU;
            inputLL(:, tt) = amp_rgc*tLL(tt).*wLL;
            inputLD(:, tt) = amp_rgc*tLD(tt).*wLD;
            inputRR(:, tt) = amp_rgc*tRR(tt).*wRR;
            inputRU(:, tt) = amp_rgc*tRU(tt).*wRU;
            inputRL(:, tt) = amp_rgc*tRL(tt).*wRL;
            inputRD(:, tt) = amp_rgc*tRD(tt).*wRD;
            
            %add external inputs onto JR at each time step tt
            
            JR = JR + inputLR(:, tt) + inputLU(:, tt) + inputLL(:, tt) + inputLD(:, tt) + inputRR(:, tt) + inputRU(:, tt) + inputRL(:, tt) + inputRD(:, tt);
            
            H = H + dt*(-H + JR )/tau; %model prediction of calcium activities at each time step for every dt.
            
            if (tLearn >= dtData) %model updates weights if tLearn exceeds dtData. Since dtData = 2*dt, this happens every other time step.
                tLearn = 0; 
                err = JR - Adata(:, iLearn+1); %As in Andalman. Adata has entries every 0.5 s and JR at every 0.25 s.
                meanerr2 = mean(err.^2); %what is displayed as chi2 when training. Use it to assess model convergence.
                chi2(nRun) = chi2(nRun) + meanerr2; %nRun is epoch number so accumulates meanerr2 every update to weight matrix. Want this to decreaase over training
                lastchi2s(iLearn) = meanerr2; %The last meanerr2 for the last weight update during an epoch.
                if ((learn) && (nRun <= nRunTot - nFree) && (nRun > nFreePre)) %nFree and nFreePre are weight freezing parameters. learn is set to 1 at start of code.
                    AR = [R(:, tt); epoch_LR; epoch_LU; epoch_LL; epoch_LD; epoch_RR; epoch_RU; epoch_RL; epoch_RD]; %augmented Dyn variable. Each trainable external input is added here.
                    %For the next 4 lines you are computing an estimate of the 
                    %inverse cross correlation matrix of network rates that
                    %are used to scale the extent of the weight update.
                    %Further details in Sussillo & Abbott (2009). 
                    k = PJ*AR; 
                    rPr = AR'*k;
                    c = 1.0/(1.0 + rPr);
                    PJ = PJ - c*(k*k');
                    %Updating external input weights if they are on
                    if epoch_LR
                        wLR = wLR - c*err*k(end-7);
                    end
                    if epoch_LU
                        wLU = wLU - c*err*k(end-6);
                    end
                    if epoch_LL
                        wLL = wLL - c*err*k(end-5);
                    end
                    if epoch_LD
                        wLD = wLD - c*err*k(end-4);
                    end
                    if epoch_RR
                        wRR = wRR - c*err*k(end-3);
                    end
                    if epoch_RU
                        wRU = wRU - c*err*k(end-2);
                    end
                    if epoch_RL
                        wRL = wRL - c*err*k(end-1);
                    end
                    if epoch_RD
                        wRD = wRD - c*err*k(end);
                    end
                    
                    J = J - c*err*k(1:N)'; %update J by err and proportional to inverse cross correlation network rates
                    
                end
                iLearn = iLearn + 1; %Set index of Adata to time of next frame
                [epoch_LR,epoch_LU,epoch_LL,epoch_LD,epoch_RR,epoch_RU,epoch_RL,epoch_RD]=deal(0); %Set epochs of external inputs to 0
            end
        end
        
        %Summary of model fit - pVar means percentage of variance explained
        rModelSample = R(:, iModelSample);
        pVar = 1 - (norm(Adata - rModelSample, 'fro')/(sqrt(N*length(tData))*stdData)).^2;
        pVars(nRun) = pVar;
        fprintf('trial=%d pVar=%f chi2=%f\n', nRun, pVar, chi2(nRun));
        
        if nRun == nRunTot
            if plotStatus
            f = figure('Position',[100 100 1800 600]);
            clf(f);
            idx = randi(N);
            figure(1);
            %subplot(5,10,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40], 'align');
            %subplot(10,10,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70], 'align');
            subplot(3,6,[1 2 3 4 5 6 7 8 9 10 11 12], 'align');
            %subplot(10,1,[1 2 3 4 5 6 7 8], 'align');
            %subplot(4,8, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24], 'align');
            hold all;
            %plot(t,R(idx,:));
            %plot(tData,Adata(idx,:));
            plot(tData(1:5:end),Adata(idx,1:5:end),'DisplayName','data');
            plot(t(1:5:end),R(idx,1:5:end),'DisplayName','model','color','r');
            lgd = legend;
            lgd.NumColumns = 1;
            lgd.Location = 'northeast';
            %lgd.Location = 'best';
            title('Example Neuron')
            txt = ['Epoch = ' int2str(nRun) ', Neuron = ' int2str(idx)];
            subtitle(txt);
            xlabel('time (s)');
            ylabel('\Deltaf/f');
            set(gca,'Box','off','TickDir','out','FontSize',9);
            %subplot(5,10,[41 42 43 44 45 46 47 48 49 50], 'align');
            subplot(3,6,[13 14 15 16 17 18], 'align');
            %subplot(10,10,[71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100], 'align');
            %subplot(10,1,[9 10], 'align');
            %subplot(4,8,[25 26 27 28 29 30 31 32], 'align');
            %hold off;
            h = heatmap(tStim(:,1:2:end));
            grid off
            axs = struct(gca); %ignore warning that this should be avoided
            cb = axs.Colorbar;
            cb.Ticks = 0:4;
            cb.TickLabels = {'No stim','Right','Forward','Left','Backward'};
            XLabels = 1:1000;
            % Convert each number in the array into a string
            CustomXLabels = string(XLabels);
            % Replace all but the fifth elements by spaces
            CustomXLabels(mod(XLabels,100) ~= 0) = " ";
            % Set the 'XDisplayLabels' property of the heatmap 
            % object 'h' to the custom x-axis tick labels
            h.XDisplayLabels = CustomXLabels;
            %YLabels = 1:2;
            %CustomYLabels(1) = "L";
            %CustomYLabels(2) = "R";
            %h.YDisplayLabels = CustomYLabels;
            h.YData = {'Left Eye','Right Eye'};
            %cmap = jet(5);
            %colormap(cmap);
            %colorbar;

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
            end

            if plotStatus
            f2 = figure('Position',[100 100 1800 600]);
            clf(f2);
            figure(2);
            %subplot(2,4,5);
            %hold off;
            plot(pVars(1:nRun));
            ylabel('pVar');
            xlabel('epoch');
            title('Percent Variance Explained');
            set(gca,'Box','off','TickDir','out','FontSize',14);
            drawnow;
            saveas(f2,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            end
            
            if plotStatus
            f3 = figure('Position',[100 100 1800 600]);
            clf(f3);
            figure(3);
            %subplot(2,4,6);
            %hold off;
            plot(chi2(1:nRun))
            ylabel('mean squared error');
            xlabel('epoch');
            title('Convergence');
            set(gca,'Box','off','TickDir','out','FontSize',14);
            drawnow;
            saveas(f3,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            end
            
            %{
            if nRun == 1
                if plotStatus
                f4 = figure('Position',[100 100 1800 600]);
                clf(f4);
                figure(4);
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
                %CustomYLabels(1) = "L";
                %CustomYLabels(2) = "R";
                Jhm.YDisplayLabels = CustomYLabels;
                %Jhm.YData = {'Left Eye','Right Eye'};
                %cmap = jet(5);
                %colormap(cmap);
                %colorbar;

                %mymap = [1 1 1
                %    1 0 0
                %    0 1 0
                %    0 0 1
                %    1 0 1];
                colormap(cmap);
                colorbar;
                c = max(abs([min(J(:,:)),max(J(:,:))]));
                caxis([-c c]);
                axs.Fontsize = 10;
                Jhm.Title = 'Synaptic Strength Weights';
                %addXLabel(Jhm,'Presynaptic','FontSize',12);
                %addYLabel(Jhm,'Postsynaptic','FontSize',12);
                drawnow;
                end
                %saveas(f4,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %saveas(f4,fullfile(fname,sprintf('epoch%d.eps',nRun)),'epsc');
                saveas(f4,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %exportgraphics(f4,fullfile(fname,sprintf('epoch%d.pdf',nRun)),'BackgroundColor','none','ContentType','vector');

            end
            %}
            %{
            if nRun == 50
                if plotStatus
                f5 = figure('Position',[100 100 1800 600]);
                clf(f5);
                figure(5);
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
                CustomXLabels(XLabels == 239) = "Bottom";
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
                %CustomYLabels(1) = "L";
                %CustomYLabels(2) = "R";
                Jhm.YDisplayLabels = CustomYLabels;
                %Jhm.YData = {'Left Eye','Right Eye'};
                %cmap = jet(5);
                %colormap(cmap);
                %colorbar;

                %mymap = [1 1 1
                %    1 0 0
                %    0 1 0
                %    0 0 1
                %    1 0 1];
                colormap(cmap);
                colorbar;
                c = max(abs([min(J(:,:)),max(J(:,:))]));
                caxis([-c c]);
                axs.Fontsize = 10;
                Jhm.Title = 'Synaptic Strength Weights';
                %addXLabel(Jhm,'Presynaptic','FontSize',12);
                %addYLabel(Jhm,'Postsynaptic','FontSize',12);
                drawnow;
                end
                %saveas(f5,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %saveas(f5,fullfile(fname,sprintf('epoch%d.eps',nRun)),'epsc');
                saveas(f5,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %exportgraphics(f5,fullfile(fname,sprintf('epoch%d.jpeg',nRun)),'BackgroundColor','none','ContentType','vector');
            end
            %}
            %{
            if nRun == 50
                if plotStatus
                f6 = figure('Position',[100 100 1800 600]);
                clf(f6);
                figure(6);
                Jhm = heatmap(J);
                S = struct(Jhm); % Undocumented
                ax = S.Axes;    % Undocumented
                % Remove grids
                Jhm.GridVisible = 'off';
                % Place lines around selected columns and row
                % Assumes columns and rows are 1 unit in size!
                col = [54 118 205 268 350];    
                row = [54 118 202 268 350];
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
                CustomXLabels(XLabels == 86) = "Right (Right Eye)";
                CustomXLabels(XLabels == 162) = "Right (Left Eye)";
                CustomXLabels(XLabels == 237) = "Backward";
                CustomXLabels(XLabels == 309) = "Left (Right Eye)";
                CustomXLabels(XLabels == 380) = "Left (Left Eye)";
                % Set the 'XDisplayLabels' property of the heatmap 
                % object 'h' to the custom x-axis tick labels
                Jhm.XDisplayLabels = CustomXLabels;
                YLabels = 1:410;
                CustomYLabels = string(YLabels);
                CustomYLabels(mod(YLabels,500) ~= 0) = " ";
                CustomYLabels(YLabels == 27) = "Forward";
                CustomYLabels(YLabels == 86) = "Right (Right Eye)";
                CustomYLabels(YLabels == 162) = "Right (Left Eye)";
                CustomYLabels(YLabels == 237) = "Backward";
                CustomYLabels(YLabels == 309) = "Left (Right Eye)";
                CustomYLabels(YLabels == 380) = "Left (Left Eye)";
                %CustomYLabels(1) = "L";
                %CustomYLabels(2) = "R";
                Jhm.YDisplayLabels = CustomYLabels;
                %Jhm.YData = {'Left Eye','Right Eye'};
                %cmap = jet(5);
                %colormap(cmap);
                %colorbar;

                %mymap = [1 1 1
                %    1 0 0
                %    0 1 0
                %    0 0 1
                %    1 0 1];
                colormap(cmap);
                colorbar;
                c = max(abs([min(J(:,:)),max(J(:,:))]));
                caxis([-c c]);
                axs.Fontsize = 10;
                Jhm.Title = 'Synaptic Strength Weights';
                %addXLabel(Jhm,'Presynaptic','FontSize',12);
                %addYLabel(Jhm,'Postsynaptic','FontSize',12);
                drawnow;
                end
                %saveas(f5,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %saveas(f5,fullfile(fname,sprintf('epoch%d.eps',nRun)),'epsc');
                saveas(f6,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
                %exportgraphics(f5,fullfile(fname,sprintf('epoch%d.jpeg',nRun)),'BackgroundColor','none','ContentType','vector');
            end
            %}
        end
        if nRun == nRunTot
            if plotStatus
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
            %CustomYLabels(1) = "L";
            %CustomYLabels(2) = "R";
            Jhm.YDisplayLabels = CustomYLabels;
            %Jhm.YData = {'Left Eye','Right Eye'};
            %cmap = jet(5);
            %colormap(cmap);
            %colorbar;

            %mymap = [1 1 1
            %    1 0 0
            %    0 1 0
            %    0 0 1
            %    1 0 1];
            colormap(cmap);
            colorbar;
            %c = max(abs([min(J(:,:)),max(J(:,:))]));
            %caxis([-c c]);
            caxis([-.65 .65]);
            axs.Fontsize = 10;
            Jhm.Title = 'Synaptic Strength Weights';
            %addXLabel(Jhm,'Presynaptic','FontSize',12);
            %addYLabel(Jhm,'Postsynaptic','FontSize',12);
            drawnow;
            saveas(f7,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            end
            %saveas(f5,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            %saveas(f5,fullfile(fname,sprintf('epoch%d.eps',nRun)),'epsc');
            %exportgraphics(f5,fullfile(fname,sprintf('epoch%d.jpeg',nRun)),'BackgroundColor','none','ContentType','vector');
        end
        %{
        if nRun == nRunTot
            if plotStatus
            f8 = figure('Position',[100 100 1800 600]);
            clf(f8);
            figure(8);
            Jhm = heatmap(J);
            S = struct(Jhm); % Undocumented
            ax = S.Axes;    % Undocumented
            % Remove grids
            Jhm.GridVisible = 'off';
            % Place lines around selected columns and row
            % Assumes columns and rows are 1 unit in size!
            col = [54 118 205 268 350];    
            row = [54 118 202 268 350];
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
            CustomXLabels(XLabels == 86) = "Right (Right Eye)";
            CustomXLabels(XLabels == 162) = "Right (Left Eye)";
            CustomXLabels(XLabels == 237) = "Backward";
            CustomXLabels(XLabels == 309) = "Left (Right Eye)";
            CustomXLabels(XLabels == 380) = "Left (Left Eye)";
            % Set the 'XDisplayLabels' property of the heatmap 
            % object 'h' to the custom x-axis tick labels
            Jhm.XDisplayLabels = CustomXLabels;
            YLabels = 1:410;
            CustomYLabels = string(YLabels);
            CustomYLabels(mod(YLabels,500) ~= 0) = " ";
            CustomYLabels(YLabels == 27) = "Forward";
            CustomYLabels(YLabels == 86) = "Right (Right Eye)";
            CustomYLabels(YLabels == 162) = "Right (Left Eye)";
            CustomYLabels(YLabels == 237) = "Backward";
            CustomYLabels(YLabels == 309) = "Left (Right Eye)";
            CustomYLabels(YLabels == 380) = "Left (Left Eye)";
            %CustomYLabels(1) = "L";
            %CustomYLabels(2) = "R";
            Jhm.YDisplayLabels = CustomYLabels;
            %Jhm.YData = {'Left Eye','Right Eye'};
            %cmap = jet(5);
            %colormap(cmap);
            %colorbar;

            %mymap = [1 1 1
            %    1 0 0
            %    0 1 0
            %    0 0 1
            %    1 0 1];
            colormap(cmap);
            colorbar;
            %c = max(abs([min(J(:,:)),max(J(:,:))]));
            %caxis([-c c]);
            caxis([-.65 .65]);
            axs.Fontsize = 10;
            Jhm.Title = 'Synaptic Strength Weights';
            %addXLabel(Jhm,'Presynaptic','FontSize',12);
            %addYLabel(Jhm,'Postsynaptic','FontSize',12);
            drawnow;
            end
            %saveas(f5,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            %saveas(f5,fullfile(fname,sprintf('epoch%d.eps',nRun)),'epsc');
            saveas(f8,fullfile(fname,sprintf('epoch%d.jpg',nRun)));
            %exportgraphics(f5,fullfile(fname,sprintf('epoch%d.jpeg',nRun)),'BackgroundColor','none','ContentType','vector');
        end
        %}
    end
    
    varData = var(reshape(Adata, N*length(tData_ind), 1));
    chi2 = chi2/(sqrt(N*length(tData_ind))*varData);
    lastchi2s = lastchi2s/(sqrt(N*length(tData_ind))*varData);
    
    %Variables you are saving from the model. You can load these variables from the
    %output in matlab and save them as a different file type.
    save(outputname, 'R', 'N', 't', 'chi2', 'Adata', 'tData', 'tData_ind', 'nRunTot', 'nFree', 'nFreePre', 'data_start_time', 'data_end_time', 'inputN', 'inputLR', 'inputLU', 'inputLL', 'inputLD', 'inputRR', 'inputRU', 'inputRL', 'inputRD', 'J', 'pVars', 'lastchi2s', 'simfilename', 'wLR', 'wLU', 'wLL', 'wLD', 'wRR', 'wRU', 'wRL', 'wRD', '-v7.3');
    
end