%%% dale_transform_example.m
   
cmap = redblue(100);
N = 100;

J = normrnd(0,1,[N,N]);
figure(1);
imshow(J);
title('Randomly initialized matrix J')

J_prime_0 = dale_transform(J, 0);
J_prime_1 = dale_transform(J, 1);

sum_of_distances = 0;

figure(2);
subplot(2, 3, 1)
h1 = heatmap(J_prime_0);
h1.Colormap = cmap;
%h1.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
caxis([-2, 2]);
h1.Title = 'RNN weight matrix J^\prime (sum-imputation)';
h1.XLabel = 'Presynaptic neuron';
h1.YLabel = 'Postsynaptic neuron';
subplot(2, 3, 2)
h2 = heatmap(J);
h2.Colormap = cmap;
%h2.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
caxis([-2, 2]);
h2.Title = 'RNN weight matrix J';
h2.XLabel = 'Presynaptic neuron';
subplot(2, 3, 3)
h3 = heatmap(J_prime_1);
h3.Colormap = cmap;
%h3.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
caxis([-2, 2]);
h3.Title = 'RNN weight matrix J^\prime (majority-imputation)';
h3.XLabel = 'Presynaptic neuron';    
subplot(2, 3, 4)
imshow(J_prime_0 == 0)
title('Zeroed Entries in J^\prime (sum-imputation)')
xlabel(sprintf('Number: %i', length(find(J_prime_0 == 0))))
subplot(2, 3, 5)
imshow(J == 0)
title('Zeroed Entries in J')
subplot(2, 3, 6)
imshow(J_prime_1 == 0)
title('Zeroed Entries in J^\prime (majority-imputation)')
xlabel(sprintf('Number: %i', length(find(J_prime_1 == 0))))

distance_0 = euclidean_distance(J, J_prime_0);
distance_1 = euclidean_distance(J, J_prime_1);
disp('distance using sum-imputation')
disp(distance_0)
disp('distance using majority-imputation')
disp(distance_1)

disp(J)
disp('Number of 0s in J')
disp(sum(sum(J == 0)))
disp('Number of 0s using sum-imputation')
disp(sum(sum(J_prime_0 == 0)))
disp('Number of 0s using majority-imputation')
disp(sum(sum(J_prime_1 == 0)))