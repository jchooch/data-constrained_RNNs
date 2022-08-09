% simple_train_holdout.m

function  [J, R, N, T, varData, chi2, pVars, val_chi2, val_pVars] = simple_train_holdout(data, inputs, max_epochs, val_epochs, lambda)
    % "NORMALIZE" DATA
    data = data/max(max(data));
    data = min(data, 0.999);
    data = max(data, -0.999);
    % SPLIT TRAIN/TEST
    train_data = data(:, 1:floor(0.8*size(data, 2)));
    val_data = data(:, floor(0.8*size(data, 2))+1:end);
    
    N = size(data, 1);
    T = size(data, 2);
    train_T = size(train_data, 2);
    val_T = size(val_data, 2);
    C = size(inputs, 1); % number of input channels
    if size(inputs, 2) ~= T
        disp(['Second dimension of inputs does not match second dimension of data. \n' ...
            'Data should have shape N X T. Inputs should have shape C X T.'])
    end
    g = 1.25;
    ampIn = 0.005;
    P0 = 1;
    amp_rgc = 0.05;
        %{
    MAYBE ADD NOISE, MAYBE NOT
    ampWN = sqrt(1/dt);
    iWN = ampWN*randn(N, length(t));
    inputN = ones(N, length(t));
    for tt = 2:length(t)
        inputN(:, tt) = iWN(:, tt) + (inputN(:, tt - 1) - iWN(:, tt))*exp(-dt); %white noise
    end
    inputN = ampIn*inputN;
    %}
    J0 = g * randn(N,N) / sqrt(N);
    J = J0;
    R = nan(N, train_T);
    AR = nan(N+8, train_T); %augmented activity for training input wts
    neuronal_inputs = zeros(C, N, train_T);
    input_weights = randn(C, N) / sqrt(N);

    PJ = P0 * eye(N+8, N+8);
    
    chi2 = zeros(max_epochs, 1);
    pVars = zeros(max_epochs, 1);
    val_chi2 = zeros(max_epochs, 1);
    val_pVars = zeros(max_epochs, 1);

    H = nan(N, T);
    
    val_counter = 0;
    for epoch = 1:max_epochs
        H(:, 1) = data(:, 1); %Initialize activities of all neurons with real values at first time point
        R(:, 1) = tanh(H(:, 1)); %nonlinearly transformed activities
        H(:, 2) = J * R(:, 1);
        for timepoint = 2:train_T-1 %time steps for each epoch. 
            R(:, timepoint) = tanh(H(:, timepoint)); %nonlinear transformation of activities
            JR = J * R(:, timepoint); % + inputN(:, timepoint); %product of synaptic strength weight matrix and nonlinear activities + gaussian noise
            for channel = 1:C
                neuronal_inputs(channel, :, timepoint) = amp_rgc * inputs(channel, timepoint) .* input_weights(channel);
                JR = JR + transpose(neuronal_inputs(channel, :, timepoint));
            end
            H(:, timepoint+1) = JR;
            err = H(:, timepoint+1) - data(:, timepoint+1) + lambda * sum(sum(abs(J))); %As in Andalman. Adata has entries every 0.5 s and JR at every 0.25 s.
            meanerr2 = mean(err.^2); %what is displayed as chi2 when training. Use it to assess model convergence.
            chi2(epoch) = chi2(epoch) + meanerr2;
            AR = [R(:, timepoint); inputs(:, timepoint)]; %augmented Dyn variable. Each trainable external input is added here.
            k = PJ*AR;
            rPr = AR'*k;
            c = 1.0/(1.0 + rPr);
            PJ = PJ - c*(k*k');
            %Updating external input weights if they are on
            for channel = 1:C
                if inputs(channel, timepoint) == 1
                        input_weights(channel, :) = input_weights(channel, :) - transpose(c * err * k(end-8-channel));
                end
            end                    
            J = J - c*err*k(1:N)'; %update J by err and proportional to inverse cross correlation network rates
        end

        R(:, T) = tanh(H(:, T));
        std_dev_data = std(reshape(data.',1,[]));
        pVar = 1 - (norm(data - R, 'fro')/(sqrt(N*T)*std_dev_data)).^2;
        pVars(epoch) = pVar;
        fprintf('trial=%d pVar=%f chi2=%f\n', epoch, pVar, chi2(epoch));

        val_neuronal_inputs = zeros(C, N, val_T);
        val_R = nan(N, val_T);
        val_H = nan(N, val_T);
        val_H(:, 1) = val_data(:, 1);
        val_R(:, 1) = tanh(val_H(:, 1));
        val_H(:, 2) = J * val_R(:, 1);
    
        for timepoint = 1:val_T-1
            real_timepoint = timepoint + train_T;
            val_R(:, timepoint) = tanh(val_H(:, timepoint));
            val_JR = J * val_R(:, timepoint);
            for channel = 1:C
                val_neuronal_inputs(channel, :, timepoint) = amp_rgc * inputs(channel, real_timepoint) .* input_weights(channel);
                val_JR = val_JR + transpose(val_neuronal_inputs(channel, :, timepoint));
            end
            val_H(:, timepoint+1) = val_JR;
            val_err = val_H(:, timepoint+1) - val_data(:, timepoint+1);
            val_meanerr2 = mean(val_err.^2);
            val_chi2(epoch) = val_chi2(epoch) + val_meanerr2;
        end
        val_R(:, val_T) = tanh(val_H(:, val_T));
        val_std_dev_data = std(reshape(val_data.',1,[]));
        val_pVar = 1 - (norm(val_data - val_R, 'fro')/(sqrt(N*val_T)*val_std_dev_data)).^2;
        val_pVars(epoch) = val_pVar;

        fprintf('trial=%d val_pVar=%f val_chi2=%f\n', epoch, val_pVar, val_chi2(epoch));

        if epoch > 1 && not(val_chi2(epoch) < val_chi2(epoch - 1))
            val_counter = val_counter + 1;
        end
        if val_counter > val_epochs
            fprintf('Validation loss stopped decreasing for %d epochs in a row, so training was stopped.\n', val_epochs)
            break
        end
    end
    varData = var(reshape(data.',1,[]));
    chi2 = chi2/(sqrt(N*T)*varData);
end