% Alexander Hay
% Homework 2

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 2\n');

%% Problem 1
% ***********************************************
% 1a

fprintf('\n');
fprintf('Problem 1a *****************************\n');
fprintf('\n');

% emulate experiment
N = 15;
mu = 6.3;
sigma = 1.9;
Ca = sigma * randn(N,1) + mu;

% standard error is sigma/root(N-1)
std_error = sigma/sqrt(N-1);

fprintf('Standard error of the mean: %.2f\n', std_error);

% ***********************************************
% 1b

fprintf('\n');
fprintf('Problem 1b *****************************\n');
fprintf('\n');

% 95% interval is mean +/- 1.96*standard error
interval_lower = mu - (1.96 * std_error);
interval_upper = mu + (1.96 * std_error);

fprintf('Lower bound: %.2f\n', interval_lower);
fprintf('Upper bound: %.2f\n', interval_upper);

% ***********************************************
% 1c

fprintf('\n');
fprintf('Problem 1c *****************************\n');
fprintf('\n');

% what is the probability that a mean of 9 is contained within the
% confidence interval
prob_1c = normcdf(9,mu,sigma);

fprintf('Probability that the true mean is 9: %.2f\n', prob_1c);
% ***********************************************
% 1d

fprintf('\n');
fprintf('Problem 1d *****************************\n');
fprintf('\n');

fprintf('Assumptions:\n');
fprintf('* Normal distribution\n');
fprintf('* Distribution is symmetric (given, if above is true)\n');
fprintf('Any variance is equal\n');

%% Problem 2
% ***********************************************
% 2a

fprintf('\n');
fprintf('Problem 2a *****************************\n');
fprintf('\n');

load('SynapseStengths_HW2.mat');

x = pulseFreq;
y = synapseStrength;

figure_2a = figure;
scatter(x,y);
title('2a - Frequency v Strength');
xlabel('Pulse Frequency (Hz)');
ylabel('Synapse Strength (pA)');

fprintf('See Figure 2a\n');

% ***********************************************
% 2b

fprintf('\n');
fprintf('Problem 2b *****************************\n');
fprintf('\n');

% arrays for part 2c
SS_res = zeros(6,1);
SS_tot = zeros(6,1);
r2 = zeros(6,1);
r2_adj = zeros(6,1);

figure_2b = figure;
scatter(x,y); %,'DisplayName','data');
title('2b - Polynomial Fitting');
xlabel('Pulse Frequency (Hz)');
ylabel('Synapse Strength (pA)');
hold on;

% fit/plot data to polynomial
% I know this isn't subplotted, I ran into a "subplot cannot convert
% handle" error :(
for i = 1:6;
    p = polyfit(x,y,i);
    y1 = polyval(p,x);
    plot(x(1:8),y1(1:8));
    
    % calculations for part c
    SS_res(i) = sum((y - y1).^2);                           % residual sum of squares
    SS_tot(i) = sum((y - mean(y)).^2 );                     % total sum of squares
    r2(i) = 1 - (SS_res(i)/SS_tot(i));                      % standard rsquared 
    r2_adj(i) = r2(i) * ((length(y)-1)/(length(y)-i));      % adjusted rsquared
end

legend_2b=cell(7,1);
legend_2b{1} = 'data';
legend_2b{2} = '1st degree';
legend_2b{3} = '2nd degree';
legend_2b{4} = '3rd degree';
legend_2b{5} = '4th degree';
legend_2b{6} = '5th degree';
legend_2b{7} = '6th degree';

legend(legend_2b,'Location','northwest');

hold off;

fprintf('See Figure 2b\n');

% ***********************************************
% 2c

fprintf('\n');
fprintf('Problem 2c *****************************\n');
fprintf('\n');

fprintf('1st degree polynomial:\n');
fprintf('R2: %.2f',r2(1)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(1));
fprintf('\n');
fprintf('2nd degree polynomial:\n');
fprintf('R2: %.2f',r2(2)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(2));
fprintf('\n');
fprintf('3rd degree polynomial:\n');
fprintf('R2: %.2f',r2(3)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(3));
fprintf('\n');
fprintf('4th degree polynomial:\n');
fprintf('R2: %.2f',r2(4)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(4));
fprintf('\n');
fprintf('5th degree polynomial:\n');
fprintf('R2: %.2f',r2(5)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(5));
fprintf('\n');
fprintf('6th degree polynomial:\n');
fprintf('R2: %.2f',r2(6)); 
fprintf('\tR2_adj: %.2f\n', r2_adj(6));

% ***********************************************
% 2d

fprintf('\n');
fprintf('Problem 2d *****************************\n');
fprintf('\n');

figure_2d = figure;
hold on;
plot(r2,'DisplayName','R^{2}');
plot(r2_adj,'DisplayName','R^{2}_{adj}');
title('2d - R^{2} & R^{2}_{adj} plots');
xlabel('Polynomial Order');
ylabel('Value');
legend();

fprintf('See Figure 2d\n');

% ***********************************************
% 2e

fprintf('\n');
fprintf('Problem 2e *****************************\n');
fprintf('\n');

fprintf('The 6th order polynomial had the highest R2 value.\n');
fprintf('The most appropriate R2 value is from the 3rd order polynomial\n');
fprintf('a lower order has a significant increase in R value\n');
fprintf('while a higher order doesnt decrease it very much\n');

%% Problem 3
% ***********************************************
% 3a

fprintf('\n');
fprintf('Problem 3a *****************************\n');
fprintf('\n');

sugar_lo = [13,19,32,34,49,15,20,19,30,8];
sugar_hi = [6,9,20,31,41,14,21,16,22,7];

fprintf('Null hypothesis: sugar does not affect the firing rates of neurons.\n');
fprintf('μ_lo = μ_hi\n');
fprintf('\n');
fprintf('Alternative hypothesis: sugar does affect the firing rate of neurons.\n');
fprintf('μ_lo ≠ μ_hi\n');

% ***********************************************
% 3b

fprintf('\n');
fprintf('Problem 3b *****************************\n');
fprintf('\n');

[h,p] = ttest2(sugar_lo,sugar_hi);

% I'm concerned I have this backwards, or straight wrong
fprintf('Test: Unpaired t-test\n');
fprintf('P value: %.3f\n',p);

% ***********************************************
% 3c

fprintf('\n');
fprintf('Problem 3c *****************************\n');
fprintf('\n');

fprintf('I think the language of the test needs to be more precise,\n');
fprintf('or the test be conducted in a broader sense, addressing\n');
fprintf('metabolism, previous diet, long term affect of high/low sugar\n');
fprintf('diets, time of sampling, etc.\n');

% ***********************************************
% 3d

fprintf('\n');
fprintf('Problem 3d *****************************\n');
fprintf('\n');

fprintf('A paired t-test should be conducted.\n');

% ***********************************************
% 3e 

fprintf('\n');
fprintf('Problem 3e *****************************\n');
fprintf('\n');

[h,p] = ttest(sugar_lo,sugar_hi);

fprintf('Test: Paired t-test\n');
fprintf('P value: %.3f\n',p);
fprintf('\n');
fprintf('I think emphasis of dietary effects bleeding over into the\n');
fprintf('next day should be more prominent. Measures should be taken\n');
fprintf('to address this.\n');

%% Problem 4
% ***********************************************
% 4a

fprintf('\n');
fprintf('Problem 4a *****************************\n');
fprintf('\n');

fprintf('Null hypothesis: protein A is not selective for inhibitory neurons.\n');
fprintf('μ_inhib = μ_excite\n');
fprintf('\n');
fprintf('Alt hypothesis: protein A is selective for inhibitatory neurons.\n');
fprintf('μ_lo ≠ μ_hi\n');
fprintf('\n');
fprintf('Assumes that if the protein does not inhibit, it excites.\n');
fprintf('The alternative would be that the protein does nothing.\n');

% ***********************************************
% 4b

fprintf('\n');
fprintf('Problem 4b *****************************\n');
fprintf('\n');

four_b = table([128;62],[35;85],'VariableNames',{'Stains for A','Does not stain for A'},'RowNames',{'Inhibitory','Excitatory'});
[h,p,stats] = fishertest(four_b);

four_b
fprintf('Fishers exact test\n');
fprintf('P value: %f\n',p);