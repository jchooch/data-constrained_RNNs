%%% mahalanobis_distance.m

function d_mahal = mahalanobis_distance(vec, vec_or_mat)
    % vec should be Nx1, vec_or_mat should be Nx1 or NxM
    vec_1 = vec;
    if size(vec_or_mat, 2) == 1
        vec_2 = vec_or_mat;
    else
        vec_2 = mean(vec_or_mat, 2);
    end
    covariance_matrix = cov(vec_or_mat.');
    inverse_covariance_matrix = inv(covariance_matrix);
    d_mahal = sqrt((vec_1 - vec_2).' * inverse_covariance_matrix * (vec_1 - vec_2));
end