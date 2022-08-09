% simple_train.m

function  [J, R, N, T, varData, chi2, pVars] = simple_train(data, inputs, number_of_epochs)
    N = size(data, 1);
    % N = double(N);
    T = size(data, 2);
    C = size(inputs, 1); % C is the number of input channels
    if size(inputs, 2) ~= T
        disp(['Second dimension of inputs does not match second dimension of data. \n' ...
            'Data should have shape N X T. Inputs should have shape C X T.'])
    end
    g = 1.25;   % Andalman used between 1.2 and 1.3 for g. Generally recommended to be between 1 and 1.5
    ampIn = 0.005; % scale of external inputs to network (h_0) (CS)
    P0 = 1; % free parameter to the learning rule for J (CS) - Recommended to be 1 to 0.01 from Sussillo & Abbott
    amp_rgc = 0.05;
    % "NORMALIZE" DATA
    data = data/max(max(data));
    data = min(data, 0.999);
    data = max(data, -0.999);
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
    R = nan(N, T);
    AR = nan(N+8, T); %augmented activity for training input wts
    neuronal_inputs = zeros(C, N, T);
    input_weights = randn(C, N) / sqrt(N);

    PJ = P0 * eye(N+8, N+8);
    
    chi2 = zeros(number_of_epochs, 1);
    pVars = zeros(number_of_epochs, 1);

    H = nan(N, T);

    for epoch = 1:number_of_epochs %Number of epochs
        H(:, 1) = data(:, 1); %Initialize activities of all neurons with real values at first time point
        R(:, 1) = tanh(H(:, 1)); %nonlinearly transformed activities
        H(:, 2) = J * R(:, 1);
        stimuli_epochs = zeros(C, 1);
        for timepoint = 2:T-1 %time steps for each epoch. 
            R(:, timepoint) = tanh(H(:, timepoint)); %nonlinear transformation of activities
            for channel = 1:C
                if inputs(channel, timepoint) == 1 
                    stimuli_epochs(channel) = 1; 
                end
            end
            JR = J * R(:, timepoint); % + inputN(:, timepoint); %product of synaptic strength weight matrix and nonlinear activities + gaussian noise
            for channel = 1:C
                neuronal_inputs(channel, :, timepoint) = amp_rgc * inputs(channel, timepoint) .* input_weights(channel);
                JR = JR + transpose(neuronal_inputs(channel, :, timepoint));
            end
            H(:, timepoint+1) = JR;
            err = JR - data(:, timepoint+1); %As in Andalman. Adata has entries every 0.5 s and JR at every 0.25 s.
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
            stimuli_epochs = zeros(C, 1);
        end
        R(:, T) = tanh(H(:, T));
        std_dev_data = std(reshape(data.',1,[]));
        pVar = 1 - (norm(data - R, 'fro')/(sqrt(N*T)*std_dev_data)).^2;
        pVars(epoch) = pVar;
        fprintf('trial=%d pVar=%f chi2=%f\n', epoch, pVar, chi2(epoch));
    end
    varData = var(reshape(data.',1,[]));
    chi2 = chi2/(sqrt(N*T)*varData);
end