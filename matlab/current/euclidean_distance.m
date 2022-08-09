% euclidean_distance.m

function dist = euclidean_distance(matrix_1, matrix_2)
    sum_of_squares = 0;
    for i = 1:length(matrix_1)
        for j = 1:length(matrix_1)
            squared_difference = (matrix_1(i,j) - matrix_2(i,j))^2;
            sum_of_squares = sum_of_squares + squared_difference;
        end
    end
    dist = sqrt(sum_of_squares);
end