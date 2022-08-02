% dale_clustering.m

number_of_nets = 100;

for filenumber = 1:number_of_nets
    load(sprintf('ensemble/output_%d.mat', filenumber), 'J');
    J = dale_transform(J, 0);
    if filenumber == 1
        accumulated_Js = J(:)';
    else
        accumulated_Js = [accumulated_Js; J(:)'];
    end
    file = filenumber;
end

fprintf('finished with filenumber %d \n', file)
fprintf('size of accumulated Js: %d \n', size(accumulated_Js))

cv = nan([1,length(J)]); % one coefficient of variation for each weight
for weight = 1:length(J)
    cv(weight) = std(J(:,weight)) / mean(J(:,weight));
end
mean_cv = sum(cv) / length(J);
fprintf('mean coefficient of variation: %d \n', mean_cv)
disp('(relative variation over networks, averaged over weights)')

% J SOLUTIONS ARE ALL OVER THE PLACE. WE CAN SEE THIS FROM:

% From t-sne with multiple distance metrics:
figure(1)
idx = kmeans(accumulated_Js, 10);
embeddings = tsne(accumulated_Js);
gscatter(embeddings(:,1), embeddings(:,2), idx)

% From t-sne with multiple distance metrics
figure(2)
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','cosine');
subplot(2,4,1);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Cosine')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','chebychev');
subplot(2,4,2);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Chebychev')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','euclidean');
subplot(2,4,3);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Euclidean')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','seuclidean');
subplot(2,4,4);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Standardised Euclidean')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','cityblock');
subplot(2,4,5);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('City Block')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','correlation');
subplot(2,4,6);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Correlation')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','spearman');
subplot(2,4,7);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Spearman')
embeddings = tsne(accumulated_Js,'Algorithm','exact','Distance','hamming');
subplot(2,4,8);
gscatter(embeddings(:,1),embeddings(:,2), idx)
title('Hamming')

% From pca
coeff = pca(accumulated_Js);
mapcaplot(accumulated_Js)

%%% Biclustering!

bicluster_result = biclust(J, 'cc'); % applies Cheng & Church algorithm 