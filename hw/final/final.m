% Alexander Hay
% Homework 6

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 6\n');

% α σ μ ≠ 

%% Problem Set 1
% ***********************************************
% 1a

fprintf('Problem Set 1\n');
fprintf('\n');

fprintf('1a\n');
fprintf('Having too many parameters to fit the data may reintroduce noise into the data youre trying to filter.\n');
fprintf('Too many parameters may also make it difficult to visualize any data clusters.\n');
fprintf('\n');

% ***********************************************
% 1b

% ***********************************************
% 1c

% ***********************************************
% 1d

fprintf('1d\n');
fprintf('Standard deviation measures the variation of the data.\n');
fprintf('Standard error measures how far the sample mean is from the true mean.\n');
fprintf('\n');

% ***********************************************
% 1e

fprintf('1e\n');
fprintf('Assumption 1a: Data is a normal distribution\n');
fprintf('Counter scenario: flu disproportionally affects children and elderly, creating two population means\n');
fprintf('\n');
fprintf('Assumption 1b: Reasonably large population (there enough data to plot a normal distribution)\n');
fprintf('Counter scenario: Sample size of one from a population of one\n');
fprintf('Test: Jarque-Bera test\n');
fprintf('\n');
fprintf('Assumption 2: Populations tested have the same variance\n');
fprintf('Counter scenario: Men and women have same average intelligence, men have more vairance in intelligence than women\n');
fprintf('Test: 2-sample F test\n');
fprintf('\n');
fprintf('Assumption 3: Data sampled is representative of the total population\n');
fprintf('Counter scenario: Testing quality of grain where the sample tested was a pebble that got mixed in with the grain\n');
fprintf('Test: Z-test\n')
fprintf('\n');

% ***********************************************
% 1f

fprintf('1f\n');
fprintf('95% confidence interval means that there is 95 certainty that the range of values contain the true mean\n');
fprintf('\n');
% ***********************************************
% 1g

fprintf('1g\n');
fprintf('p-values indicate the strength of the evidence against the null-hypothesis. A p-value of 0.05 shows a 5% probability\n');
fprintf('that the null hypothesis is correct\n');
fprintf('\n');

%% Problem Set 2
% ***********************************************
% 2a

rat = linspace(1,15,15);
lo_anx = [3,6,3,7,7,2,3,5,8,10,2,5,7,8,11];
hi_anx = [5,6,7,11,8,4,8,9,12,17,15,3,16,14,9];

fprintf('2a ********** \n');
fprintf('Low anxiety\n');
fprintf('mean:\t%.2f \t sigma:\t%.2f\n',mean(lo_anx),std(lo_anx));
fprintf('\n');
fprintf('High anxiety\n');
fprintf('mean:\t%.2f \t sigma:\t%.2f\n',mean(hi_anx),std(hi_anx));
fprintf('\n');

%% **********************************************
% 2b

n = length(rat);
mean_lo_anx = mean(lo_anx);
mean_hi_anx = mean(hi_anx);
std_lo_anx = std(lo_anx);
std_hi_anx = std(hi_anx);

% standard error is sigma/root(N-1) (per hw2)
sem_lo_anx = std_lo_anx/sqrt(n-1);
sem_hi_anx = std_hi_anx/sqrt(n-1);

lo_anx_int_lo = mean_lo_anx - (1.96 * sem_lo_anx);
lo_anx_int_hi = mean_lo_anx + (1.96 * sem_hi_anx);

hi_anx_int_lo = mean_hi_anx - (1.96 * sem_lo_anx);
hi_anx_int_hi = mean_hi_anx + (1.96 * sem_hi_anx);

fprintf('2b ********** \n');
fprintf('Low anxiety\n');
fprintf('95%% interval: %.2f - %.2f\n', lo_anx_int_lo, lo_anx_int_hi);
fprintf('\n');
fprintf('High anxiety\n');
fprintf('95%% interval: %.2f - %.2f\n', hi_anx_int_lo, hi_anx_int_hi);
fprintf('\n');

%% **********************************************
% 2c

fprintf('2c ********** \n');
fprintf('A two sample t-test would determine if the data comes from different populations or not.\n');
fprintf('In this context it would determine if each condition have different means (effects).\n');
fprintf('\n');

%% **********************************************
% 2d

fprintf('2d ********** \n');
fprintf('Null hypothesis: that both conditions have the same means\n');
fprintf('μ_lo_anx = μ_hi_anx\n')
fprintf('\n');
fprintf('Alt hypothesis: that both conditions have the different means\n');
fprintf('μ_lo_anx ≠ μ_hi_anx\n')
fprintf('\n');

%% **********************************************
% 2e

[h p ci stats] = ttest2(lo_anx, hi_anx);

fprintf('2e ********** \n');
fprintf('Two-sample t-test results:\n');
fprintf('h:\t%.f\n',h);
fprintf('DoF:\t%.f\n',stats.df);
fprintf('p:\t%.4f\n',p);
fprintf('\n');

%% **********************************************
% 2f

fprintf('2f ********** \n');
fprintf('The results of the two-sample t-test indicate that the data from each condition represent different effects\n');
fprintf('\n');

%% **********************************************
% 2g

fprintf('2g ********** \n');
fprintf('Statistical power:\t\t%.2f\n',sampsizepwr('t2',[mean_lo_anx std_lo_anx],mean_hi_anx,[],n));
fprintf('Rats needed for power of 0.95:\t%.f\n',sampsizepwr('t2',[mean_lo_anx std_lo_anx],mean_hi_anx,.95,[]));
fprintf('Statistical power is the likelihood of the test rejecting the hypothesis\n');
fprintf('A power of 0.95 demonstrates a 95%% confidence of rejecting the null hypothesis\n');
fprintf('\n');

%% **********************************************
% 2h

fprintf('2h ********** \n');
fprintf('Using the same rats would mean using a paired t-test\n');
fprintf('Statistical power:\t%.2f\n',sampsizepwr('t',[mean_lo_anx std_lo_anx],mean_hi_anx,[],n));
fprintf('Advantages:\trequires fewer rats\n');
fprintf('\t\teliminates variations between rats\n');
