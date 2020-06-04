% Alexander Hay
% Homework 3

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 3\n');

% α σ μ ≠ 
%% Problem 1
% ***********************************************
% 1a

fprintf('\n');
fprintf('Problem 1a *****************************\n');
fprintf('\n');

n = 5;                              % number of rats
spd = [4,5,7,2,2];                  % seizures per day

mu = sum(spd)/n
sigma = sqrt((spd-mu)*(spd-mu)'/n);
sigma_fitdist = sqrt((spd-mu)*(spd-mu)'/(n-1));

pd = fitdist(spd','Normal');

fprintf('Calculated mean: \t%.2f\n', mu);
fprintf('Calculated std dev: \t%.2f\n', sigma);
fprintf('\n');
fprintf('fitdist mean: \t\t%.2f\n', pd.mu);
fprintf('fitdist std dev: \t%.2f\n',pd.sigma);
fprintf('\n');
fprintf('They have different values, because fistdist uses maximum likelyhood\n');
fprintf('See code for calculation difference.\n')
fprintf('\n');
fprintf('σ = \t\t%.2f\n', sigma)
fprintf('σ_fitdist = \t%.2f\n', sigma_fitdist);

% ***********************************************
% 1b

fprintf('\n');
fprintf('Problem 1b *****************************\n');
fprintf('\n');

phat = mle(spd);
fprintf('MLE mean: \t%.2f\n', phat(1));
fprintf('MLE std dev: \t%.2f\n', phat(2));
fprintf('\n');
fprintf('They have the same values, because MLE calculates sigma the same way\n');
fprintf('\n');
fprintf('σ = \t\t%.2f\n', sigma);
fprintf('σ_MLE = \t%.2f\n', phat(2));

% Use parameter values for HW
sigma = phat(2);

% ***********************************************
% 1c

fprintf('\n');
fprintf('Problem 1c *****************************\n');
fprintf('\n');

n_new = sampsizepwr('t', [mu sigma], 2, .90, [], 'tail', 'left');

fprintf('Samples needed: %.0f\n', n_new);

% ***********************************************
% 1d

fprintf('\n');
fprintf('Problem 1d *****************************\n');
fprintf('\n');

n_new = sampsizepwr('t', [mu sigma], 2, .99, [], 'tail', 'left');

fprintf('Samples needed: %.0f\n', n_new);

% ***********************************************
% 1e

fprintf('\n');
fprintf('Problem 1e *****************************\n');
fprintf('\n');

mu_new = sampsizepwr('t', [mu sigma], [], .90, 6, 'tail', 'left');

fprintf('New mean: %.2f\n', mu_new);


% ***********************************************
% 1f

fprintf('\n');
fprintf('Problem 1f *****************************\n');
fprintf('\n');

p_1 = sampsizepwr('t', [mu sigma], 1, [], 6, 'tail', 'left');
p_2 = sampsizepwr('t', [mu sigma], 3, [], 6, 'tail', 'left');

fprintf('Statistical power (1/day): %.2f\n', p_1);
fprintf('Statistical power (3/day): %.2f\n', p_2);

% ***********************************************
% 1g

fprintf('\n');
fprintf('Problem 1g *****************************\n');
fprintf('\n');

n_new = sampsizepwr('t2', [mu sigma], 2, .9, [], 'tail', 'left');

fprintf('Samples needed: %.0f\n', n_new);

%% Problem 2
% ***********************************************
% 2a

fprintf('\n');
fprintf('Problem 2a *****************************\n');
fprintf('\n');

load('SynapseData.mat');

% synapseData = synapseData';         % puts data in better form         

% (trace location, time)

%    | t1 | t2 | t3 | ...
%    +----+----+----+-----
% l1 | mV | mV | mV | ...
% l2 | mV | mV | mV | ...
% l3 | mV | mV | mV | ...
% ...| ...| ...| ...| ...

trace1 = synapseData(1,:);

figure_2a = figure;

% label = 'Trace: ' + string(t);
plot(trace1,'DisplayName','Trace 1');
hold on;

title('2a - First trace location');
xlabel('time (ms)');
xlim([0 500]);
xticklabels({'0','10','20','30','40','50'});
ylabel('trace voltage (mV)');
legend();

alpha = max(trace1); % amplitude
j = 50;      % how the line drops
k = 250;   % decay ("flexiness")

% make array for V(t)
V_size = size(synapseData);
% V = zeros(V_size(1),1);
% V = zeros(V_size(1),V_size(2));

% fill array

for t = 1:V_size(2)
    V(t) = alpha * ((1/(1+exp(-t)))+exp((-t-j)/k)-1);
end

test_SS_res = sum((trace1 - V).^2);                           % residual sum of squares
test_SS_tot = sum((trace1 - mean(trace1)).^2 );                     % total sum of squares
test_r2 = 1 - (test_SS_res/test_SS_tot);
test_r2_adj = test_r2 * ((length(trace1)-1)/(length(trace1)-100)); 

plot(V,'DisplayName','V(t)');
hold off;

fprintf('see figure 2a\n');
fprintf('\n');
fprintf('α: %.1f\n',alpha);
fprintf('d: %.0f\n',j);
fprintf('τ: %.0f\n',k);

% ***********************************************
% 2b

fprintf('\n');
fprintf('Problem 2b *****************************\n');
fprintf('\n');

fprintf('there is no 2b');
fprintf('\n');

% ***********************************************
% 2c

fprintf('\n');
fprintf('Problem 2c *****************************\n');
fprintf('\n');

fprintf('Takes about 5s to run\n');

% create arrays to hold data
[n_trace t_samps] = size(synapseData);
r = zeros(n_trace,20,50);
r_max = zeros(100,4);

% these loops and functions essentially brute force a solution
% α is pinned to max amplitude value, to decrease runtime
% d is sampled between 0 and 200, incremented by 10
% τ is sampled between 0 and 500, incremented by 10

% each trace has its v calculated, then r_adjusted calculated
% once each variable has been sampled, the largest r value is found
% with corresponding d and τ values. The information is stored in a
% 100x4 matrix:

%        | α | d | τ | r |
%        +---+---+---+---|
% trace1 | . | . | . | . |
% trace2 | . | . | . | . |
% trace3 | . | . | . | . |
%   ...  |...|...|...|...|

for n = 1:n_trace
    trace = synapseData(n,:);
    alpha = max(trace);
    for j = 1:20
        for k = 1:50
            for t = 1:t_samps
                v(t) = v_calc(t,alpha,j,k);
            end
            r(n,j,k) = r_calc(trace,v);
        end
    end
%     fprintf('%.0f\n',n);
    r_max(n,1) = alpha;
    [r_max(n,4),d,tau] = find_max(squeeze(r(n,:,:)));
    r_max(n,2) = d*10;
    r_max(n,3) = tau*10;
end

r_max_2c = r_max;

% fprintf('Why did I spend so much time on this?\n');
fprintf('See code for details. Information is stored in r_max_2c\n');

% ***********************************************
% 2d

fprintf('\n');
fprintf('Problem 2d *****************************\n');
fprintf('\n');

mean_alpha = sum(r_max(:,1))/n_trace;
mean_d = sum(r_max(:,2))/n_trace;
mean_tau = sum(r_max(:,3))/n_trace;

std_alpha = sqrt((r_max(:,1)-mean_alpha)'*(r_max(:,1)-mean_alpha)/n_trace);
std_d = sqrt((r_max(:,2)-mean_d)'*(r_max(:,2)-mean_d)/n_trace);
std_tau = sqrt((r_max(:,3)-mean_tau)'*(r_max(:,3)-mean_tau)/n_trace);

% outliers defined as 1.5 standard deviations away from the mean
for n = 1:n_trace
    temp_alpha = r_max(n,1);
    temp_d = r_max(n,2);
    temp_tau = r_max(n,3);
    if temp_alpha < (mean_alpha - std_alpha)
        % flag for removal
        r_max(n,:) = 0;
    elseif temp_alpha > (mean_alpha + std_alpha)
        r_max(n,:) = 0;
    elseif temp_d < (mean_d - std_d)    
        r_max(n,:) = 0;
    elseif temp_d > (mean_d + std_d)
        r_max(n,:) = 0;
    elseif temp_tau < (mean_tau - std_tau)    
        r_max(n,:) = 0;
    elseif temp_tau > (mean_tau + std_tau)
        r_max(n,:) = 0;
    end
end

fprintf('Outliers defined as 1.5 standard deviations away from the mean.\n');
fprintf('They were found using a series of if/then statements.\n');

% ***********************************************
% 2e

fprintf('\n');
fprintf('Problem 2e *****************************\n');
fprintf('\n');

% remove zeros from array
r_max(r_max==0)=[];

% zeros were put in deliberately in the previous step
% because if the index was removed mid loop, the array
% may become shorter than the iteration (index) of the loop

% reshape is then used to convert back into a 2D array
% divided by 4 because r_max has/had 4 columns[α,d,τ,r]

r_max = reshape(r_max,[length(r_max)/4,4]);

h_alpha = chi2gof(r_max(:,1));
h_d = chi2gof(r_max(:,2));
h_tau = chi2gof(r_max(:,3));

figure2e = figure;

subplot2e = subplot(1,3,1,'Parent',figure2e);
histogram(r_max(:,1));
title('α histogram');
xlabel('α values');
ylabel('# of instances');
ylim([0 30]);

subplot2e = subplot(1,3,2,'Parent',figure2e);
histogram(r_max(:,2));
title('d histogram');
xlabel('d values');
ylim([0 30]);

subplot2e = subplot(1,3,3,'Parent',figure2e);
histogram(r_max(:,3));
title('τ histogram');
xlabel('τ values');
ylim([0 30]);

fprintf('Chi-square goodness-of-fit\n');
fprintf('h_α: \t%.0f\n', h_alpha);
fprintf('h_d: \t%.0f\n', h_d);
fprintf('h_τ: \t%.0f\n', h_tau);
fprintf('\n');
fprintf('Histograms do not approach a normal distribution.\n');
fprintf('α does because it was pinned to max amplitude and outliers removed.\n');
fprintf('(ie. the data was trimmed to fit a normal distribution)\n');

%% Problem 3
% ***********************************************
% 3a

fprintf('\n');
fprintf('Problem 3a *****************************\n');
fprintf('\n');

load('HW3_FMRI_data.mat');

voxels_mean = mean(voxels_faceImage);

figure3a = figure;
scatter(voxels_baseline,voxels_mean,'DisplayName','Sample');
hold on;
plot([0:1200],'DisplayName','1:1 reference line');
title('3a - face voxel baseline vs mean');
xlabel('Voxel baseline');
ylabel('Voxel mean');
xlim([-50 1200]);
ylim([-50 1200]);
legend('Location','northwest');
hold off;

fprintf('see figure 3a\n');
fprintf('\n');
fprintf('Voxel response initially appears distributed.\n');
fprintf('Most, if not all, appear to respond to the image.\n');
fprintf('But this depends what defines a response (ie. magnitude).\n');

% ***********************************************
% 3b

fprintf('\n');
fprintf('Problem 3b *****************************\n');
fprintf('\n');

voxels_test = voxels_faceImage(1:3,:);

% h_test_1 = chi2gof(voxels_test(1,:));
% h_test_2 = chi2gof(voxels_test(2,:));
% h_test_3 = chi2gof(voxels_test(3,:));

tt2_12 = ttest2(voxels_test(1,:),voxels_test(2,:));
tt2_13 = ttest2(voxels_test(1,:),voxels_test(3,:));
tt2_23 = ttest2(voxels_test(2,:),voxels_test(3,:));

fprintf('Two-sample t-test for samples 1 & 2: %.0f\n',tt2_12);
fprintf('Two-sample t-test for samples 1 & 3: %.0f\n',tt2_13);
fprintf('Two-sample t-test for samples 2 & 3: %.0f\n',tt2_23);

fprintf('Samples tested reject the hypothesis that the samples came from\n');
fprintf('distributions with different means.\n');

% ***********************************************
% 3c

fprintf('\n');
fprintf('Problem 3c *****************************\n');
fprintf('\n');

fprintf('Im going to use the fitdist function to fit the data to a normal.\n');
fprintf('distribution, then a ttest with α=0.05. The two-sample t-test\n');
fprintf('results provide a good foundation for this test.\n');

% ***********************************************
% 3d

fprintf('\n');
fprintf('Problem 3d *****************************\n');
fprintf('\n');

voxels_gaus = fitdist(voxels_mean','Normal');

% voxels_ttest = zeros(50,2);

for i = 1:1000
    [h p] = ttest(voxels_faceImage(:,i),voxels_gaus.mu,'Alpha',0.05);
    voxels_ttest(i,1) = h;
    voxels_ttest(i,2) = p;
end

figure3d = figure;

semilogy(voxels_ttest(:,2));
hold on;
semiline = [0:1000, 0.05];
plot(semiline);
title('3d - p value semilogy α=0.05');
xlabel('voxel');
ylabel('p-value');
xlim([0 1000]);
hold off;

fprintf('see figure 3d\n');

% ***********************************************
% 3e

fprintf('\n');
fprintf('Problem 3e *****************************\n');
fprintf('\n');

fprintf('Not at all. I dont think my math is off but I think I used the wrong test.\n');

% ***********************************************
% 3f

fprintf('\n');
fprintf('Problem 3f *****************************\n');
fprintf('\n');

for i = 1:1000
    [h p] = ttest(voxels_faceImage(:,i),voxels_gaus.mu,'Alpha',(0.05/1000));
    voxels_ttest(i,1) = h;
    voxels_ttest(i,2) = p;
end

figure3f = figure;

semilogy(voxels_ttest(:,2));
hold on;
semiline = [0:1000, (0.05/1000)];
plot(semiline);
title('3f - p value semilogy α=0.05/1000');
xlabel('voxel');
ylabel('p-value');
xlim([0 1000]);
hold off;

fprintf('see figure 3f\n');
fprintf('Same results as 3d, strengthening the suspicion that I used the wrong test.\n');

% ***********************************************
% 3g

fprintf('\n');
fprintf('Problem 3g *****************************\n');
fprintf('\n');

fprintf('N/A\n');

fprintf('\n');
fprintf('Notes *****************************\n');
fprintf('\n');
fprintf('I understand the format of this assignment. It was walking through\n');
fprintf('how we would use these tools IRL. I got stuck on 2c for awhile, more\n');
fprintf('of a coding issue than conceptual, but I worked through it. However,\n');
fprintf('once I started working through question 3 I felt like I needed more\n');
fprintf('...something. So far this class has been a lot on how to do the math\n');
fprintf('and this was the first problem to test why we do the math. This would\n');
fprintf('be nice to work through in class.\n');