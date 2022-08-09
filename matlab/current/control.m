% control.m

g = 1.25;   % Andalman used between 1.2 and 1.3 for g. Generally recommended to be between 1 and 1.5
dtData = 0.5; % ? frame rate given in Supplemental Info (CS) - 1/frame rate - our actual frame rate is slightly faster than 2 so actual dtData is somewhat lower
dt = 0.25; % integration step in sec - this value was used in Andalman
ampIn = 0.005; % scale of external inputs to network (h_0) (CS)
P0 = 1; % free parameter to the learning rule for J (CS) - Recommended to be 1 to 0.01 from Sussillo & Abbott

cmap = redblue(100);

N = 100; % 100 neurons
C = 8; % 8 input channels
J_dgRNN = normrnd(0,1,[N,N]); % target weight matrix
saved_input = matfile('inputs.mat');
inputs = saved_input.inputs;
T = size(inputs, 2);
I = binornd(1, 0.5, C, T); % random input

figure(1)
imshow(I(:,1:100), 'InitialMagnification', 2000)
title('Toy binary inputs (first 100/5000 timepoints)')
xlabel('Time')
ylabel('Input channel')
axis on
set(gca,'box','off')

figure(2)
h = heatmap(J_dgRNN);
h.Colormap = cmap;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h.Title = 'RNN weight matrix J_{target}';
h.XLabel = 'Presynaptic neuron';
h.YLabel = 'Postsynaptic neuron';

[H, R_targets, N_output] = dgRNN_fwd(J_dgRNN, T, I);

figure(3)
imshow(R_targets, 'InitialMagnification', 2000)
title('Neuronal activations during forward pass')
xlabel('Time')
ylabel('Neuron')

R_bio = load('calact410Nstd.mat');
R_bio = R_bio.calact410Nstd;
%[J_dcRNN, R_model, N, T, varData, chi2, pVars] = dcRNN_train(R_targets, inputs, 50, 0);
[J_dcRNN, R_model, N, T, varData, chi2, pVars] = dcRNN_train(R_bio, inputs, 50, 0);

figure(4)
plot(chi2)

figure(5)
plot(pVars)

figure(6)
subplot(3, 1, 1)
title('Neuron 1')
plot(1:1:100, R_model(1, 1:100), 1:1:100, R_targets(1, 1:100))
legend('model', 'target')
subplot(3, 1, 2)
title('Neuron 50')
plot(500:1:600, R_model(50, 500:600), 500:1:600, R_targets(50, 500:600))
legend('model', 'target')
subplot(3, 1, 3)
title('Neuron 100')
plot(1000:1:1100, R_model(100, 1000:1100), 1000:1:1100, R_targets(100, 1000:1100))
legend('model', 'target')

%{
for filenumber = 1:100
    joes_main_rnn_code(filenumber, false)
end
%}