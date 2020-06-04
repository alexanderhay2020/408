% Alexander Hay
% Homework 5

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 5\n');
fprintf('\n');
fprintf('submitted late due to network card failure in my laptop\n');


% α σ μ ≠ 
%% Problem 1
% ***********************************************
% 1a

fprintf('\n');
fprintf('Problem 1a *****************************\n');
fprintf('\n');

load('hw5_pca.mat');
% rows - trials (10,000)
% cols - samples (1,000)

figure_1a = figure;
subplot(1,2,1,'Parent',figure_1a);
imagesc(D);
% cmap = [0 0 1; 1 1 1; 1 0 0];
% cmap_interp = interp1(-1:1:1, cmap, -1:.01:1);
% colormap(cmap_interp);
title('1a - Heatmap of D');
xlabel('time sample');
ylabel('trial');

subplot(1,2,2,'Parent',figure_1a);
imagesc(D(5770:5775,:));

colorbar();
title('1a - Subsamples of D');
xlabel('time sample');
mV = colorbar;
ylabel(mV, 'mV')
yticklabels({'5770','5771','5772','5773','5774','5775'});

fprintf('see figure 1a');
fprintf('\n');

%% ***********************************************
% 1b

fprintf('\n');
fprintf('Problem 1b *****************************\n');
fprintf('\n');

figure_1b = figure;
imagesc(cov(D));
colorbar();
title('1b - Covariance matrix of D');
xlabel('time sample');
ylabel('time sample');

fprintf('see figure 1b\n');
fprintf('\n');
fprintf('The covariance is a measure of how the voltage sample of each nerve vary together.\n');
fprintf('Positive values mean the samples positively vary together, negative values mean\n');
fprintf('the samples inversely vary together.\n');
fprintf('There seems to be a time dependence going on with the voltage.\n');

%% ***********************************************
% 1c

fprintf('\n');
fprintf('Problem 1c *****************************\n');
fprintf('\n');

[coeff,score,eigens] = pca(D);
fracVar = eigens./sum(eigens);
figure_1c = figure;
subplot(1,2,1,'Parent',figure_1c);
plot(fracVar, 'bx');
title('1c - Principal Component Analysis');
xlabel('Eigenvalue number (cropped)');
ylabel('Fraction variance explained');
xlim([0 10]); % x-axis cropped to view significant eigenvalues better

subplot(1,2,2,'Parent',figure_1c);
plot(cumsum(fracVar));
title('Cumulative Sum of Fraction of Variance');
xlabel('Eigenvalue number');
ylabel('Fraction variance explained');

fprintf('see figure 1c\n');
fprintf('\n');
fprintf('After 3 PCs there is no significant contribution to explain the data.\n');
fprintf('\n');
fprintf('Variance explained:\t %.2f%% \n',sum(fracVar(1:4))*100);

%% ***********************************************
% 1d

fprintf('\n');
fprintf('Problem 1d *****************************\n');
fprintf('\n');

figure_1d = figure;
hold on;
plot(coeff(:,1),'b');
plot(coeff(:,2),'r');
plot(coeff(:,3),'g');
% plot(coeff(:,4),'k');
title('1d - Principal Components 1-3');
xlabel('time sample');
ylabel('mV');
legend('PC1', 'PC2', 'PC3');
hold off;

fprintf('see figure 1d\n');
fprintf('\n');
fprintf('The data represents reasonable EPSPs.\n');
fprintf('The 4th PC looks very noisy. Uncomment the line to show the noise.\n');

%% ***********************************************
% 1e

fprintf('\n');
fprintf('Problem 1e *****************************\n');
fprintf('\n');

figure_1e = figure;
hold on;
histogram(score(:,1),'facecolor','b');
histogram(score(:,2),'facecolor','r');
histogram(score(:,3),'facecolor','g');
histogram(score(:,4),'facecolor','k');
title('1e - Histogram of Scores');
xlabel('sample');
ylabel('score');
xlim([-500 1000]);
legend('PC score 1','PC score 2','PC score 3','PC score 4');
hold off;

fprintf('see figure 1e\n');
fprintf('\n');
fprintf('The histogram shows that the 4th PC did not contribute significantly.\n');
fprintf('The scores did not have the same magnitude as the first three.\n');

%% ***********************************************
% 1f

fprintf('\n');
fprintf('Problem 1f *****************************\n');
fprintf('\n');

figure_1f = figure;
subplot(1,3,1,'Parent',figure_1f);
scatter(score(:,1), score(:,2));
title('1f- PC1 v PC2');
xlabel('PC1 score');
ylabel('PC2 score');
subplot(1,3,2,'Parent',figure_1f);
scatter(score(:,1), score(:,3));
title('PC1 v PC3');
xlabel('PC1 score');
ylabel('PC3 score');
subplot(1,3,3,'Parent',figure_1f);
scatter(score(:,2), score(:,3));
title('PC2 v PC32');
xlabel('PC2 score');
ylabel('PC3 score');

fprintf('see figure 1f\n');
fprintf('\n');
fprintf('There are 8 clusters\n');

%% ***********************************************
% 1g

fprintf('\n');
fprintf('Problem 1g *****************************\n');
fprintf('\n');

idx = kmeans([score(:,1),score(:,2),score(:,3)],8);

figure_1g = figure;
scatter3(score(:,1), score(:,2), score(:,3),20, idx, 'filled', 'MarkerEdgeColor', 'k');
title('1g - Kmeans Classification');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

% find(idx==x) returns indicies of cluster x
% D(find(idx==x) returns array populated by cluster x
% mean(D(find(idx==x) returns mean of that array

means = [mean(D(find(idx == 1)))
         mean(D(find(idx == 2)))
         mean(D(find(idx == 3)))
         mean(D(find(idx == 4)))
         mean(D(find(idx == 5)))
         mean(D(find(idx == 6)))
         mean(D(find(idx == 7)))
         mean(D(find(idx == 8)))];

stds = [std(D(find(idx == 1)))
        std(D(find(idx == 2)))
        std(D(find(idx == 3)))
        std(D(find(idx == 4)))
        std(D(find(idx == 5)))
        std(D(find(idx == 6)))
        std(D(find(idx == 7)))
        std(D(find(idx == 8)))];

figure_1g2 = figure;
subplot(1,2,1,'Parent',figure_1g2);
errorbar([1:8],means,stds);
hold on;
title('1g - Mean Traces of Each Cluster');
xlabel('Cluster');
ylabel('mV');
xlim([0 9]);

subplot(1,2,2,'Parent',figure_1g2);
bar([1:8],means);
title('1g - Mean Traces of Each Cluster');
xlabel('Cluster');
ylabel('mV');
xlim([0 9]);

fprintf('see figures 1g\n');

%% ***********************************************
% 1h

fprintf('\n');
fprintf('Problem 1h *****************************\n');
fprintf('\n');

fprintf('This explains the 3 PCs found in part d.\n');
fprintf('The other clusters likely represent noise or residual waveforms.\n');


%% ***********************************************
% 1i

fprintf('\n');
fprintf('Problem 1i *****************************\n');
fprintf('\n');


[temp,originalpos] = sort(abs(means), 'descend' );
n = temp(1:3);
p=originalpos(1:3);

fprintf('Event rate of each synapse:\n');
fprintf('\n');
for i=1:3
    fprintf('Synapse %i:\t%.0f\n',p(i),length(find(idx==p(i))));
end

%% Problem 2
% ***********************************************
% 2a

fprintf('\n');
fprintf('Problem 2a *****************************\n');
fprintf('\n');

load('HW5_machineLearning.mat');

r_count = 0;
l_count = 0;

% Determines how much data there is for each eye
for i=1:209
    if cell2mat(whichEye(i))=='R'
        r_count = r_count+1;
    elseif cell2mat(whichEye(i))=='L'
        l_count = l_count+1;
    end
end

% Creates arrays for each eye
r_array = zeros(r_count,1133);
l_array = zeros(l_count,1133);
r_index = 1;
l_index = 1;

% Separates data for each eye
for i=1:209
    if cell2mat(whichEye(i))=='R'
        r_array(r_index,:) = transcriptomics_data(i,:);
        r_index = r_index + 1;
    elseif cell2mat(whichEye(i))=='L'
        l_array(l_index,:) = transcriptomics_data(i,:);
        l_index = l_index + 1;
    end
end

fprintf('Cells for right eye: \t %.0f\n',r_count);
fprintf('Cells for left eye: \t%.0f\n',l_count);

%% ***********************************************
% 2b

fprintf('\n');
fprintf('Problem 2b *****************************\n');
fprintf('\n');

r_mean = mean(r_array);
l_mean = mean(l_array);

exp_ratio = r_mean./l_mean;

fprintf('see workspace variables r_mean, l_mean, and exp_ratio.\n');

%% ***********************************************
% 2c

fprintf('\n');
fprintf('Problem 2c *****************************\n');
fprintf('\n');

figure_2c = figure;
histogram(exp_ratio);
title('2c - Expression Ratio');
xlabel('Ratio (mean(R)/mean(L)');
ylabel('Value');
xlim([0 2.2]);

fprintf('Histogram values around 1 indicate similar gene expressions.\n');
fprintf('Lower histogram values indicate a greater left gene expression, while higher\n');
fprintf('histogram values indicate greater right gene expression.\n');

%% ***********************************************
% 2d

fprintf('\n');
fprintf('Problem 2d *****************************\n');
fprintf('\n');

pd = fitdist(exp_ratio','Normal');

fprintf('Mean: \t%.2f\n',pd.mu);
fprintf('Sigma: \t%.2f\n',pd.sigma);

%% ***********************************************
% 2e

fprintf('\n');
fprintf('Problem 2e *****************************\n');
fprintf('\n');

prob_val_hi = normcdf(2,pd.mu,pd.sigma,'upper');
prob_val_lo = normcdf(0.5,pd.mu,pd.sigma);

fprintf('Probability value >2.0: \t%f\n',prob_val_hi);
fprintf('Probability value <0.5: \t%f\n',prob_val_lo);

%% ***********************************************
% 2f

fprintf('\n');
fprintf('Problem 2f *****************************\n');
fprintf('\n');

index_lo = find(exp_ratio <= 0.5);
index_hi = find(exp_ratio >= 2);

x = [index_lo index_hi];
y = [exp_ratio(index_lo) exp_ratio(index_hi)];

figure_2f = figure;
bar(x,y);
title('2f - Expression Ratio Bar Graph');
xlabel('Genes');
ylabel('Expression Ratio');

fprintf('Granted theres only 5, but the expression ratios are very high, then very low shortly after.\n');
fprintf('This suggests that they may be time dependent.\n');

%% ***********************************************
% 2g

fprintf('\n');
fprintf('Problem 2g *****************************\n');
fprintf('\n');

T = [array2table(transcriptomics_data) whichEye'];

fprintf('see variable T in the workspace.\n');
fprintf('\n');
fprintf('Variables:\t1,133\n');
fprintf('Classes:\t2\n');
%% ***********************************************
% 2h

fprintf('\n');
fprintf('Problem 2h *****************************\n');
fprintf('\n');

fprintf('Best Classifier:\tEnsemble\n');
fprintf('Accuracy:\t\t62.7%% \n');
fprintf('Error:\t\t\t%.1f%%\n',(78/209)*100);
fprintf('L every time:\t\t%.1f%%\n',(l_count/209)*100);
fprintf('\n');
fprintf('The classifier is barely better than just choosing L\n');
fprintf('\n');
fprintf('was this supposed to be with a data table made from part f?\n');
fprintf('Because it was made using the full dataset.\n');
fprintf('Either way its marginally better than choosing L all the time and flipping a coin for the difference.\n');

%% ***********************************************
% 2i

fprintf('\n');
fprintf('Problem 2i *****************************\n');
fprintf('\n');

index_lo = find(exp_ratio <= 1/1.3);
index_hi = find(exp_ratio >= 1.3);

newData = [transcriptomics_data(:,index_lo) transcriptomics_data(:,index_hi)];

T2 = [array2table(newData) whichEye'];

fprintf('Genes within 1.3 threshold:\t%.f\n',(numel(index_lo) + numel(index_hi)));
fprintf('\n');
fprintf('Classifier:\tQuadratic SVM\n');
fprintf('\n');
fprintf('Performance was increased due to a reduction in noise, or rather, only including\n');
fprintf('meaningful data.\n');

%% ***********************************************
% 2j

fprintf('\n');
fprintf('Problem 2j *****************************\n');
fprintf('\n');

fprintf('Incorrectly labeling L eye being from R eye:\t%.2f%%\n',(15/233)*100);
fprintf('Incorrectly labeling R eye being from L eye:\t%.2f%%\n',(20/233)*100);

%% ***********************************************
% 2k

fprintf('\n');
fprintf('Problem 2k *****************************\n');
fprintf('\n');

value = 1.37;
% 1.2, 79.4
% 1.25, 81.8
% 1.3, 84.7
% 1.32, 82.8
% 1.35, 85.2
% 1.36, 84.2
% 1.37, 86.6
% 1.38, 84.7
% 1.4, 83.3
% 1.5, 82.8

index_lo = find(exp_ratio <= 1/value);
index_hi = find(exp_ratio >= value);

newData = [transcriptomics_data(:,index_lo) transcriptomics_data(:,index_hi)];

T2 = [array2table(newData) whichEye'];

fprintf('I was able to get 86.6%% accuracy using a ratio of 1.37, using SVMs.\n');
fprintf('Accuracy deviated with each training with each session, indicating a random influence?\n');
fprintf('Id have to look at my notes, but values between 1.35 and 1.4 performed well.\n');