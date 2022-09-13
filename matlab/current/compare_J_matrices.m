%%% compare_J_matrices.m

function [J1_zeros, J2_zeros, J1_tiny, J2_tiny, J1_mean, J2_mean, J1_var, J2_var, P_man, P_aut] = compare_J_matrices(J1, J2, soft_sparsity_threshold)
    fprintf('Soft sparsity threshold: %f\n', soft_sparsity_threshold)

    % compare sparsity
    J1_zeros = numel(J1) - nnz(J1);
    J2_zeros = numel(J2) - nnz(J2);
    J1_tiny = sum(abs(J1) < soft_sparsity_threshold, 'all');
    J2_tiny = sum(abs(J2) < soft_sparsity_threshold, 'all');
    fprintf('Hard sparsity: J1: %d (%d %%). J2: %d (%d %%). \n', ...
        J1_zeros, (J1_zeros / numel(J1)), J2_zeros, (J2_zeros / numel(J2)))
    fprintf('Soft sparsity: J1: %d (%d %%).  J2: %d (%d %%). \n', ...
        J1_tiny, (J1_tiny / numel(J1)), J2_tiny, (J2_tiny / numel(J2)))

    % means and variances
    J1_mean = mean(J1, 'all');
    J2_mean = mean(J2, 'all');
    J1_var = var(J1, 0, 'all');
    J2_var = var(J2, 0, 'all');
    fprintf('Means: J1: %f. J2: %f.\n', J1_mean, J2_mean)
    fprintf('Variances: J1: %f. J2: %f.\n', J1_var, J2_var)

    % pointwise distances and overall distance
    % D = pointwise distance matrix
    D = abs(J1 - J2);
    D_mean = mean(D, 'all');
    fprintf('Mean pointwise difference: %d \n', D_mean)

    % look at mean absolute values
    % 
end