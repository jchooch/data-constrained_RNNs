% dgRNN_fwd.m
% An in silico toy data-generating RNN
% Inputs: weight_matrix, inputs, number_of_timesteps
% Outputs: H=raw_activities, R=tanh_activities, N=number_of_neurons

function [H, R, N] = dgRNN_fwd(J, T, I, W, dt, tau)
    N = length(J);
    if ~exist('dt', 'var')
        dt = 1;
    end
    if ~exist('tau', 'var')
        tau = 1;
    end
    if ~exist('I', 'var') % check if there are inputs (i.e. stimuli)
        disp('No inputs provided, so running without inputs.')
        with_inputs = false;
    elseif size(I, 2) ~= T % check if inputs are right size wrt T
        disp(['Second dimension of inputs is longer than number of timesteps. \n ' ...
            'Inputs should have dimensions NUMBER_OF_INPUT_CHANNELS X T.'])
        return % halt if inputs don't fit T
    else % i.e. if there are inputs with the right number of timesteps...
        with_inputs = true;
        number_of_input_channels = size(I, 1);
        if ~exist('W', 'var')
            W = normrnd(0,1,[N, number_of_input_channels]); % create random input weight matrix [ALTERNATIVELY, something like deal(randn(N,1)/sqrt(N)) ...?]
        elseif size(W, 1) ~= N || size(W, 2) ~= number_of_input_channels
            disp(['Weight matrix for external inputs (W) has wrong shape. \n' ...
                'It should have shape N X NUMBER_OF_INPUT_CHANNELS.'])
            return % halt if W is wrong size
        end
    end
    if length(size(J)) ~= 2 || size(J, 1) ~= size(J, 2)
        disp('Input weight matrix is wrong size. It must be 2D and square.')
        return   % halt if J is wrong size
    end
    weighted_inputs = W * I;
    
    %{
    % ADD SOME NOISE HERE AND UNCOMMENT NOISE ADDITION IN FORWARD PASS
    noise_mu = 0;
    noise_sigma = 0.01;
    input_noise = normrnd(noise_mu, noise_sigma, N, T);
    %}

    H = nan(N, T);
    R = nan(N, T);
    H(:, 1) = normrnd(0,1,[N, 1]); % randomly initialize neuron activities from standard normal distribution
    R(:, 1) = tanh(H(:, 1));

    if with_inputs == true
        H(:, 2) = J * R(:, 1) + weighted_inputs(:, 1); % + input_noise(:, 1);
        for timepoint = 2:T
            R(:, timepoint) = tanh(H(:, timepoint));
            H(:, timepoint+1) = J * R(:, timepoint) + weighted_inputs(:, timepoint); % + input_noise(:, timepoint);
        end 
    elseif with_inputs == false
        H(:, 2) = J * R(:, 1); % + input_noise(:, 1);
        for timepoint = 2:T
            R(:, timepoint) = tanh(H(:, timepoint));
            H(:, timepoint+1) = J * R(:, timepoint); % + input_noise(:, timepoint);
        end
    end
    fprintf('Completed forward pass with N=%d neurons over T=%d timesteps \n', N, T)
end