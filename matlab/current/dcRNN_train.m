%%% dcRNN_train.m

function  [J, R, N, T, varData, chi2, pVars] = dcRNN_train(data, inputs, number_of_epochs, error_mode, lambda)
    N = size(data, 1);
    T = size(data, 2);
    C = size(inputs, 1); % C is the number of input channels
    % size(inputs, 2) is number of model timepoints
    %{
    if size(inputs, 2) ~= T
        disp(['Second dimension of inputs does not match second dimension of data. \n' ...
            'Data should have shape N X T. Inputs should have shape C X T.'])
    end
    %}
    
    g = 1.25;   % Andalman used between 1.2 and 1.3 for g. Generally recommended to be between 1 and 1.5
    dtData = 0.5; % ? frame rate given in Supplemental Info (CS) - 1/frame rate - our actual frame rate is slightly faster than 2 so actual dtData is somewhat lower
    dt = 0.25; % integration step in sec - this value was used in Andalman
    ampIn = 0.005; % scale of external inputs to network (h_0) (CS)
    P0 = 1; % free parameter to the learning rule for J (CS) - Recommended to be 1 to 0.01 from Sussillo & Abbott
    %amp - Appears Andalman et al. recommends 10x larger than amplitude of noise
    amp_rgc = 0.05;
    
    tData = [1:dtData:length(data)*dtData]; %How many times J is updated

    t = 0:dt:1000;  %Number of integration time steps
    
    iModelSample = zeros(length(tData), 1);
    for i=1:length(tData)
        [tModelSample, iModelSample(i)] = min(abs(tData(i)-t));
    end
    
    data = data/max(max(data));
    data = min(data, 0.999);
    data = max(data, -0.999);
    
    %stdData = std(reshape(data, [N*length(tData_ind), 1]));
    
    %Gaussian white noise with amplitude no more than 1% of baseline
    %according to Andalman
    %set seed - have random number gen and then randn
    ampWN = sqrt(1/dt);
    iWN = ampWN*randn(N, length(t));
    inputN = ones(N, length(t));
    
    for tt = 2: length(t)
        inputN(:, tt) = iWN(:, tt) + (inputN(:, tt - 1) - iWN(:, tt))*exp(-dt); %white noise
    end
    inputN = ampIn*inputN;
    
    %Initialize 
    J0 = g*randn(N,N)/sqrt(N); %synaptic strength matrix
    J = J0;
    R = NaN(N, length(t)); %nonlinear activity
    AR = NaN(N+8, length(t)); %augmented activity for training input wts
    JR = zeros(N, 1); %product of synaptic strength weights and nonlinear activity

    neuronal_inputs = zeros(C, N, T);
    input_weights = randn(C, N) / sqrt(N);
    
    PJ = P0 * eye(N+8, N+8);
    
    data_steps = size(data, 2);
    fprintf('data_steps: %d\n', data_steps);
    fprintf('dt: %f\n', dt);
    fprintf('length(t): %d\n', length(t));
    assert(data_steps >= length(t) * dt);

    chi2 = zeros(1, number_of_epochs);
    pVars = zeros(1, number_of_epochs);

    H = nan(N, T);

    for epoch = 1:number_of_epochs %Number of epochs
        H(:, 1) = data(:, 1); %Initialize activities of all neurons with real values at first time point
        R(:, 1) = tanh(H(:, 1)); %nonlinearly transformed activities
        H(:, 2) = J * R(:, 1);
        tLearn = 0; %param for when to update J matrix
        iLearn = 1; %Used to index Adata to subtract from model predicted rates in err.
        
        stimuli_epochs = zeros(C, 1);

        for tt = 2:length(t)-2 %time steps for each epoch. 
            tLearn = tLearn + dt; %update for each time step
            R(:, tt) = tanh(H(:, tt)); %nonlinear transformation of activities
            for channel = 1:C
                if inputs(channel, tt) == 1 
                    stimuli_epochs(channel) = 1; 
                end
            end
            JR = J * R(:, tt) + inputN(:, tt); %product of synaptic strength weight matrix and nonlinear activities + gaussian noise
            % size is still good here
            %external inputs at each time step - input on if time vectors are 1 at time tt
            for channel = 1:C
                neuronal_inputs(channel, :, tt) = amp_rgc * inputs(channel, tt) .* input_weights(channel);
            end

            for channel = 1:C
                JR = JR + transpose(neuronal_inputs(channel, :, tt));
            end
            
            % H(:, tt+1) = H(:, tt) + dt*(-H(:, tt) + JR); %model prediction of calcium activities at each time step for every dt.
            % H(:, tt+1) = JR;
            H(:, tt+1) = H(:, tt) + dt * (-H(:, tt) + JR);

            if (tLearn >= dtData) %model updates weights if tLearn exceeds dtData. Since dtData = 2*dt, this happens every other time step.
                tLearn = 0;
                if error_mode == 1 % non-squared difference, no L1 penalty
                    err = JR - data(:, iLearn+1); %As in Andalman. Adata has entries every 0.5 s and JR at every 0.25 s.
                elseif error_mode == 2 % non-squared difference, L1 penalty: all weights
                    err = JR - data(:, iLearn+1) + lambda * sum(abs(J), 'all');
                elseif error_mode == 3 % non-squared difference, L1 penalty: postsynaptic weights
                    err = JR - data(:, iLearn+1) + lambda * sum(abs(J)).';
                elseif error_mode == 4 % non-squared difference, sign-sensitive L1 penalty: all weights                        
                    err = JR - data(:, iLearn+1);
                    signs = sign(err);
                    penalties = signs * lambda * sum(abs(J), 'all');
                    err = err + penalties;
                elseif error_mode == 5 % squared difference, no L1 penalty
                    err = (JR - data(:, iLearn+1)).^2;
                elseif error_mode == 6 % squared difference, L1 penalty: all weights
                    err = (JR - data(:, iLearn+1)).^2 + lambda * sum(abs(J), 'all');
                elseif error_mode == 7 % squared difference, L1 penalty: postsynaptic weights
                    err = (JR - data(:, iLearn+1)).^2 + lambda * sum(abs(J)).';
                end
                meanerr2 = mean(err.^2); %what is displayed as chi2 when training. Use it to assess model convergence.
                chi2(epoch) = chi2(epoch) + meanerr2;
                if (epoch <= number_of_epochs) && (epoch > 0)
                    AR = [R(:, tt); stimuli_epochs]; %augmented Dyn variable. Each trainable external input is added here.
                    % Estimate inverse cross correlation matrix of network rates that are used to scale
                    % the weight update. Further details in Sussillo & Abbott (2009).
                    k = PJ*AR; 
                    rPr = AR'*k;
                    c = 1.0/(1.0 + rPr);
                    PJ = PJ - c*(k*k');
                    %Updating external input weights if they are on
                    for channel = 1:C
                        if stimuli_epochs(channel) == 1
                            input_weights(channel, :) = input_weights(channel, :) - transpose(c * err * k(end-8-channel));
                        end
                    end                    
                    J = J - c*err*k(1:N)'; %update J by err and proportional to inverse cross correlation network rates
                end
                iLearn = iLearn + 1; %Set index of Adata to time of next frame
                stimuli_epochs = zeros(C, 1);
            end
        end
        
        %Summary of model fit - pVar means percentage of variance explained
        rModelSample = R(:, iModelSample);
        %pVar = 1 - (norm(data - rModelSample, 'fro')/(sqrt(N*length(tData))*stdData)).^2;
        std_dev_data = std(reshape(data.',1,[]));
        pVar = 1 - (norm(data - R(:, 1:2:end-1), 'fro')/(sqrt(N*length(tData))*std_dev_data)).^2; % used to be just R
        pVars(epoch) = pVar;
        fprintf('trial=%d pVar=%f chi2=%f\n', epoch, pVar, chi2(epoch));
    end
    
    varData = var(reshape(data, N*T, 1));
    chi2 = chi2/(sqrt(N*T)*varData);
end