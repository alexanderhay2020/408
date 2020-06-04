%% load data
load('spikeWaveformData');

%% plot example waveform
figure(1);
timeAxis = 0:5E-2:4;
plot(timeAxis, spikeWaveforms(1,:));
xlabel('Time (ms)');
ylabel('Voltage (uV)');

%% plot heatmap of all waveforms
figure(2);
imagesc(spikeWaveforms);
cmap = [0 0 1; 1 1 1; 1 0 0];
cmap_interp = interp1(-1:1:1, cmap, -1:.01:1);
colormap(cmap_interp);
caxis([-600 600]);
c = colorbar;
c.Label.String = 'Voltage';
xticks(1:10:81);
xticklabels(timeAxis(1:10:81));
xlabel('Time (ms)');
ylabel('Cell number');

%% label each cell with its type
yticks(1:10:400);
yticklabels(types_2classes(1:10:400));
ylabel('Cell type');

%% Examine covariance matrices
% covarance across cells
C = cov(spikeWaveforms');
figure(3);
imagesc(C);
xlabel('Cells');
ylabel('Cells');
title('Covariance of data comparing cells');
cmap = [0 0 1; 1 1 1; 1 0 0];
cmap_interp = interp1(-1:1:1, cmap, -1:.01:1);
colormap(cmap_interp);
caxis([-3000 3000]);
c = colorbar;
c.Label.String = 'Covariance';

%% covariance across time
C = cov(spikeWaveforms);
figure(4);
imagesc(C);
xlabel('time sample');
ylabel('time sample');
title('Covariance of data comparing time points');
cmap = [0 0 1; 1 1 1; 1 0 0];
cmap_interp = interp1(-1:1:1, cmap, -1:.01:1);
colormap(cmap_interp);
caxis([-2000 2000]);
c = colorbar;
c.Label.String = 'Covariance';

%% See what PCA tells us
[V, score, eigvals] = pca(spikeWaveforms);
fractionVariance = eigvals./sum(eigvals);
figure(5);
plot(fractionVariance, 'bx');
xlabel('Eigenvalue number');
ylabel('Fraction variance explained');

%% See how much variance is explained
sum(fractionVariance(1:4))

%% plot the first 4 eigenvectors
figure(6);
hold on;
plot(timeAxis, V(:,1),'b');
plot(timeAxis, V(:,2),'r');
plot(timeAxis, V(:,3),'g');
plot(timeAxis, V(:,4),'k');
xlabel('Time (ms)');
ylabel('Weight in eigenvector');
title('First 4 eigenvectors (PCs)');
legend('PC1', 'PC2', 'PC3', 'PC4');
hold off;

%% some 2D scatter plots
typeColors = strcmp(types_2classes, 'TypeA');

figure(7);
subplot(2,3,1);
scatter(score(:,1), score(:,2), [], typeColors);
xlabel('PC1 score');
ylabel('PC2 score');
subplot(2,3,2);
scatter(score(:,1), score(:,3), [], typeColors);
xlabel('PC1 score');
ylabel('PC3 score');
subplot(2,3,3);
scatter(score(:,1), score(:,4), [], typeColors);
xlabel('PC1 score');
ylabel('PC4 score');
subplot(2,3,4);
scatter(score(:,2), score(:,3), [], typeColors);
xlabel('PC2 score');
ylabel('PC3 score');
subplot(2,3,5);
scatter(score(:,2), score(:,4), [], typeColors);
xlabel('PC2 score');
ylabel('PC4 score');
subplot(2,3,6);
scatter(score(:,3), score(:,4), [], typeColors);
xlabel('PC3 score');
ylabel('PC4 score');

%% Assemble table for classification learner
T = [array2table(spikeWaveforms) types_2classes];

%% Trained the classifier from the command line
[trainedClassifier, validationAccuracy] = trainClassifier(T)

%% Load new data set
load('../mat_files/spikeWaveformData_newData');

%% Predict classes for new data
predictions = trainedClassifier.predictFcn(array2table(spikeWaveforms))

%% now try 4 classes
load('../mat_files/spikeWaveformData');
T = [array2table(spikeWaveforms) types_4classes];
typeA1 = strcmp(types_4classes, 'TypeA1');
typeA2 = strcmp(types_4classes, 'TypeA2');
typeB1 = strcmp(types_4classes, 'TypeB1');
typeB2 = strcmp(types_4classes, 'TypeB2');
typeColors = typeA1*1 + typeA2*2 + typeB1*3 + typeB2*4;

figure(8);
subplot(2,3,1);
scatter(score(:,1), score(:,2), [], typeColors);
xlabel('PC1 score');
ylabel('PC2 score');
subplot(2,3,2);
scatter(score(:,1), score(:,3), [], typeColors);
xlabel('PC1 score');
ylabel('PC3 score');
subplot(2,3,3);
scatter(score(:,1), score(:,4), [], typeColors);
xlabel('PC1 score');
ylabel('PC4 score');
subplot(2,3,4);
scatter(score(:,2), score(:,3), [], typeColors);
xlabel('PC2 score');
ylabel('PC3 score');
subplot(2,3,5);
scatter(score(:,2), score(:,4), [], typeColors);
xlabel('PC2 score');
ylabel('PC4 score');
subplot(2,3,6);
scatter(score(:,3), score(:,4), [], typeColors);
xlabel('PC3 score');
ylabel('PC4 score');

