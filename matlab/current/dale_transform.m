% dale_transform.m

% First arg is matrix to be Dale-transformed.
% Second arg is either 0 or 1, depending on how
%   E/I identity of presynaptic neurons should be imputed.
%   0: Sum-imputation takes the identity of the neuron to be the sign of
%   the sum of the weights.
%   1: Majority-imputation takes the identity of the neuron to be the
%   majority-sign of the weights.
%   So 'sum' takes into account weight magnitudes, but 'majority' doesn't.
%   Can add more identity_imputation settings: random, experimental, 
%   majority-above-threshold, etc.
function J_prime = dale_transform(J, identity_imputation)
    if identity_imputation == 0
        J_prime = J;
        identities = zeros([1, length(J)]);
        for i = 1:length(J_prime) % loop through columns
            if sum(J_prime(:,i)) >= 0
                identities(i) = 1;
            end
            if sum(J_prime(:,i)) < 0
                identities(i) = -1;
            end
            for j = 1:length(J_prime)
                if identities(i) < 0 && J_prime(j,i) >= 0
                    J_prime(j,i) = 0;
                end
                if identities(i) >= 0 && J_prime(j,i) < 0
                    J_prime(j,i) = 0;
                end
            end
        end
    end
    if identity_imputation == 1
        J_prime = J;
        identities = zeros([1, length(J)]);
        for i = 1:length(J_prime) % loop through presynaptic neurons
            pos_count = 0;
            neg_count = 0;
            for j = 1:length(J_prime) % loop through output weights of the ith presynaptic
                if J_prime(j, i) >= 0
                    pos_count = pos_count + 1; % count negative weights
                end
                if J_prime(j, i) < 0
                    neg_count = neg_count + 1; % count positive weights
                end
            end
            if pos_count >= neg_count   % set identity based on majority sign
                identities(i) = 1;
            else
                identities(i) = -1;
            end
        end
        for i = 1:length(J_prime)
            for j = 1:length(J_prime)
                if identities(i) > 0 && J_prime(j,i) < 0
                    J_prime(j,i) = 0;
                end
                if identities(i) < 0 && J_prime(j,i) >= 0
                    J_prime(j,i) = 0;
                end
            end
        end
    end
end