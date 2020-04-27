% Alexander Hay
% Homework 1

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 1\n');

%% Problem 1

A=[5,4,4,3,4,5,4,3,5,3,2,2,6,5,4,4,4,4,2,5];
B=[6,6,7,5,6,5,7,6,5,7,8,8,7,6,5,6,5,6,4,9];
C=[5,3,2,4,6,4,3,4,5,4,1,2,2,2,3,2,6,5,3,5];

% ***********************************************
% 1a

% find unique values for each gene
a=sort(unique(A));
b=sort(unique(B));
c=sort(unique(C));

fprintf('\n');
fprintf('Problem 1a *****************************\n');
fprintf('\n');
fprintf('Values in gene A, sorted:\n');
fprintf('[%i %i %i %i %i]\n', a);
fprintf('Values in gene B, sorted:\n');
fprintf('[%i %i %i %i %i %i]\n', b);
fprintf('Values in gene B, sorted:\n');
fprintf('[%i %i %i %i %i %i]\n', c);
fprintf('\n');

% ***********************************************
% 1b

% Find number of instances of each number
% Put it in array
for i = 2:6
    P_A(i-1) = sum(A(:)==i);
end

for i = 4:9
    P_B(i-3) = sum(B(:)==i);
end

for i = 1:6
    P_C(i) = sum(C(:)==i);
end

% number of instances of each number
% divided by number of samples
prob_A = P_A/numel(A);
prob_B = P_B/numel(B);
prob_C = P_C/numel(C);

fprintf('\n');
fprintf('Problem 1b *****************************\n');
fprintf('\n');
fprintf('Probability check for A\n');
fprintf('Sample means of A:\n');
fprintf('[   %i    %i    %i    %i    %i]\n',a);
fprintf('[%.2f %.2f %.2f %.2f %.2f]\n',prob_A);
fprintf('\n');
fprintf('Sample means of B:\n');
fprintf('[   %i    %i    %i    %i    %i    %i]\n',b);
fprintf('[%.2f %.2f %.2f %.2f %.2f %.2f]\n',prob_B);
fprintf('\n');
fprintf('Sample means of C:\n');
fprintf('[   %i    %i    %i    %i    %i    %i]\n',c);
fprintf('[%.2f %.2f %.2f %.2f %.2f %.2f]\n',prob_C);

% ***********************************************
% 1c

prob_A_check = sum(prob_A);
prob_B_check = sum(prob_B);
prob_C_check = sum(prob_C);

fprintf('\n');
fprintf('Problem 1c *****************************\n');
fprintf('\n');
fprintf('Probability check for A\n');
fprintf('%.2f\n', prob_A_check);
fprintf('\n');
fprintf('Probability check for B\n');
fprintf('%.2f\n', prob_B_check);
fprintf('\n');
fprintf('Probability check for C\n');
fprintf('%.2f\n', prob_C_check);

% ***********************************************
% 1d

fprintf('\n');
fprintf('Problem 1d *****************************\n');
fprintf('\n');
fprintf('see figure 1\n');

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(prob_A);
ylabel('Probability');
xlabel('Values');
title('Problem 1d');
ylim(axes1,[0 1]);
box(axes1,'on');
set(axes1,'XTick',[1 2 3 4 5],'XTickLabel',{'2','3','4','5','6'});

% ***********************************************
% 1e

mean_A = sum(A)/length(A);
mean_B = sum(B)/length(B);
mean_C = sum(C)/length(C);

fprintf('\n');
fprintf('Problem 1e *****************************\n');
fprintf('\n');
fprintf('Sample mean for A\n');
fprintf('%.2f\n', mean_A);
fprintf('\n');
fprintf('Sample mean for B\n');
fprintf('%.2f\n', mean_B);
fprintf('\n');
fprintf('Sample mean for C\n');
fprintf('%.2f\n', mean_C);
fprintf('\n');

% ***********************************************
% 1f

% variance is the average of sample-mean difference, squared, 
% for all samples
var_A = (A-mean_A)*(A-mean_A)'/length(A);
var_B = (B-mean_B)*(B-mean_B)'/length(B);
var_C = (C-mean_C)*(C-mean_C)'/length(C);

% std dev is the square root of the variance
std_dev_A = sqrt(var_A);
std_dev_B = sqrt(var_B);
std_dev_C = sqrt(var_C);

fprintf('\n');
fprintf('Problem 1f *****************************\n');
fprintf('\n');
fprintf('Variance and standard deviation for A\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_A, std_dev_A);
fprintf('\n');
fprintf('Variance and standard deviation for B\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_B, std_dev_B);
fprintf('\n');
fprintf('Variance and standard deviation for C\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_C, std_dev_C);
fprintf('\n');

%% Problem 2
% ***********************************************
% 2a

X1 = [1,2,3,4,5,6];
X2 = [1,2,3,4,5,6];
fprintf('\n');
fprintf('Problem 2a *****************************\n');
fprintf('\n');
% fprintf('X1:\n');
% fprintf('%i %i %i %i %i %i\n',X1);
% fprintf('X2:\n');
% fprintf('%i %i %i %i %i %i\n',X2);
% fprintf('\n');
for i = 1:length(X1)
    for j = 1:length(X2)
        Y_mat(i,j) = [X1(i)+X2(j)];
    end
end

% Possile values of random variable
% Number of combinations to create random variable

P = [unique(Y_mat)';
    [sum(Y_mat(:)==2),
     sum(Y_mat(:)==3),
     sum(Y_mat(:)==4),
     sum(Y_mat(:)==5),
     sum(Y_mat(:)==6),
     sum(Y_mat(:)==7),
     sum(Y_mat(:)==8),
     sum(Y_mat(:)==9),
     sum(Y_mat(:)==10),
     sum(Y_mat(:)==11),
     sum(Y_mat(:)==12)]'];

Y = P(1,:);
 
true_prob = P(2,:)/numel(Y_mat);

fprintf('[Possile values of random variable]\n');
fprintf('Y = [  %i    %i    %i    %i    %i    %i    %i    %i    %i   %i   %i]\n',Y);
fprintf('\n');
fprintf('[Number of combinations to create random variable]\n');
fprintf('P = [%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f]\n',true_prob);
fprintf('\n');

% ***********************************************
% 2b
 
% takes number of possilbe combinations to
% create variable and divides it by total number
% of possible combinations
fprintf('\n');
fprintf('Problem 2b *****************************\n');
fprintf('\n');

true_prob = P(2,:)/numel(Y_mat);
true_prob_check = sum(true_prob);

fprintf('true probabilities:\n');
fprintf('P = [%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f]\n',true_prob);
fprintf('\n');
fprintf('sum of probabilities: %f\n',true_prob_check);

% ***********************************************
% 2c

figure2 = figure;
axes2 = axes('Parent',figure2);
hold(axes1,'on');
bar(true_prob,'DisplayName','true_prob');
title({'Problem 2c'});
ylabel({'Probability'});
ylim([0 1]);
xlabel({'Possible Y Values'});
box(axes2,'on');
set(axes2,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});

fprintf('\n');
fprintf('Problem 2c *****************************\n');
fprintf('\n');
fprintf('see figure 2\n');


% ***********************************************
% 2d

fprintf('\n');
fprintf('Problem 2d *****************************\n');
fprintf('\n');

mean_Y = sum(Y)/length(Y);

fprintf('True mean of random variable Y\n');
fprintf('%.2f\n', mean_Y);

% ***********************************************
% 2e

fprintf('\n');
fprintf('Problem 2e *****************************\n');
fprintf('\n');

var_Y = (Y-mean_Y)*(Y-mean_Y)'/length(Y);
std_dev_Y = sqrt(var_Y);

fprintf('Variance and standard deviation for Y\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_Y, std_dev_Y);

% ***********************************************
% 2f

fprintf('\n');
fprintf('Problem 2f *****************************\n');
fprintf('\n');

mean_X1 = sum(X1)/length(X1);
mean_X2 = sum(X2)/length(X2);
var_X1 = (X1-mean_X1)*(X1-mean_X1)'/length(X1);
var_X2 = (X2-mean_X2)*(X2-mean_X2)'/length(X2);
std_dev_X1 = sqrt(var_X1);
std_dev_X2 = sqrt(var_X2);

fprintf('Variance and standard deviation for X1\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_X1, std_dev_X1);
fprintf('\n');
fprintf('Variance and standard deviation for X2\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_X2, std_dev_X2);
fprintf('\n');
fprintf('variance and std dev are compounded\n');

% ***********************************************
% 2g

fprintf('\n');
fprintf('Problem 2g *****************************\n');
fprintf('\n');

% Fano-factor = std dev^2/mean
FF_Y = (std_dev_Y^2)/mean_Y;

% Coefficient of variation = std dev/mean
CV_Y = std_dev_Y/mean_Y;

fprintf('Fano-factor and coefficient of variance for Y\n');
fprintf('Fano Factor: %.2f | Coefficient of Variance: %.2f\n', FF_Y, CV_Y);
fprintf('\n');

% ***********************************************
% 2h

fprintf('\n');
fprintf('Problem 2g *****************************\n');
fprintf('\n');
fprintf('Multinomial\n');

%% Problem 3
% ***********************************************
% 3a

fprintf('\n');
fprintf('Problem 3a *****************************\n');
fprintf('\n');

X1_100 = randi(6,100,1);
X2_100 = randi(6,100,1);
Y_100 = X1_100 + X2_100;

fprintf('The array is too large to print\n');
fprintf('the array is saved as trials_1000 in the workspace\n');

% ***********************************************
% 3b 

fprintf('\n');
fprintf('Problem 3b *****************************\n');
fprintf('\n');

X1_100 = randi(6,100,1);
X2_100 = randi(6,100,1);
Y_100 = X1_100 + X2_100;

% find probabilities of each value
trials_100_probs = [sum(Y_100(:)==2)/100,
                    sum(Y_100(:)==3)/100,
                    sum(Y_100(:)==4)/100,
                    sum(Y_100(:)==5)/100,
                    sum(Y_100(:)==6)/100,
                    sum(Y_100(:)==7)/100,
                    sum(Y_100(:)==8)/100,
                    sum(Y_100(:)==9)/100,
                    sum(Y_100(:)==10)/100,
                    sum(Y_100(:)==11)/100,
                    sum(Y_100(:)==12)/100];

figure3 = figure;
axes4 = axes('Parent',figure3);
hold(axes4,'on');
bar(trials_100_probs);
ylabel({'Probability'});
xlabel('Y Value (X1 + X2)');
title({'Problem 3b - 100 Trials'});
xlim(axes4,[-0.2 12]);
ylim(axes4,[0 1]);
set(axes4,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});

fprintf('see figure 3\n');


% ***********************************************
% 3c

fprintf('\n');
fprintf('Problem 3c *****************************\n');
fprintf('\n');

X1_1000 = randi(6,1000,1);
X2_1000 = randi(6,1000,1);
Y_1000 = X1_1000 + X2_1000;

% find probabilities of each value
trials_1000_probs = [sum(Y_1000(:)==2)/1000,
                     sum(Y_1000(:)==3)/1000,
                     sum(Y_1000(:)==4)/1000,
                     sum(Y_1000(:)==5)/1000,
                     sum(Y_1000(:)==6)/1000,
                     sum(Y_1000(:)==7)/1000,
                     sum(Y_1000(:)==8)/1000,
                     sum(Y_1000(:)==9)/1000,
                     sum(Y_1000(:)==10)/1000,
                     sum(Y_1000(:)==11)/1000,
                     sum(Y_1000(:)==12)/1000];

X1_10000 = randi(6,10000,1);
X2_10000 = randi(6,10000,1);
Y_10000 = X1_10000 + X2_10000;

% find probabilities of each value
trials_10000_probs = [sum(Y_10000(:)==2)/10000,
                      sum(Y_10000(:)==3)/10000,
                      sum(Y_10000(:)==4)/10000,
                      sum(Y_10000(:)==5)/10000,
                      sum(Y_10000(:)==6)/10000,
                      sum(Y_10000(:)==7)/10000,
                      sum(Y_10000(:)==8)/10000,
                      sum(Y_10000(:)==9)/10000,
                      sum(Y_10000(:)==10)/10000,
                      sum(Y_10000(:)==11)/10000,
                      sum(Y_10000(:)==12)/10000];                 

X1_100000 = randi(6,100000,1);
X2_100000 = randi(6,100000,1);
Y_100000 = X1_100000 + X2_100000;

% find probabilities of each value
trials_100000_probs = [sum(Y_100000(:)==2)/100000,
                       sum(Y_100000(:)==3)/100000,
                       sum(Y_100000(:)==4)/100000,
                       sum(Y_100000(:)==5)/100000,
                       sum(Y_100000(:)==6)/100000,
                       sum(Y_100000(:)==7)/100000,
                       sum(Y_100000(:)==8)/100000,
                       sum(Y_100000(:)==9)/100000,
                       sum(Y_100000(:)==10)/100000,
                       sum(Y_100000(:)==11)/100000,
                       sum(Y_100000(:)==12)/100000];                  
figure4 = figure;

subplot1 = subplot(2,2,1,'Parent',figure4);
hold(subplot1,'on');
bar(trials_100_probs);
title('100 trials');
box(subplot1,'on');
set(subplot1,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});
ylabel({'Probability'});
xlabel({'Value'});
ylim([0 1]);
axes4 = axes('Parent',figure4);

subplot2 = subplot(2,2,2,'Parent',figure4);
hold(subplot2,'on');
bar(trials_1000_probs);
title('1,000 trials');
box(subplot2,'on');
set(subplot2,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});
ylabel({'Probability'});
xlabel({'Value'});
ylim([0 1]);

subplot3 = subplot(2,2,3,'Parent',figure4);
hold(subplot3,'on');
bar(trials_10000_probs);
title('10,000 trials');
box(subplot3,'on');
set(subplot3,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});
ylabel({'Probability'});
xlabel({'Value'});
ylim([0 1]);

subplot4 = subplot(2,2,4,'Parent',figure4);
hold(subplot4,'on');
bar(trials_100000_probs);
title('100,000 trials');
box(subplot4,'on');
set(subplot4,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});
ylabel({'Probability'});
xlabel({'Value'});
ylim([0 1]);

fprintf('see figure 4\n');

% ***********************************************
% 3d

fprintf('\n');
fprintf('Problem 3d *****************************\n');
fprintf('\n');

mean_Y_100 = sum(Y_100)/length(Y_100);
mean_Y_1000 = sum(Y_1000)/length(Y_1000);
mean_Y_10000 = sum(Y_10000)/length(Y_10000);
mean_Y_100000 = sum(Y_100000)/length(Y_100000);

var_Y_100 = (Y_100-mean_Y_100)'*(Y_100-mean_Y_100)/length(Y_100);
var_Y_1000 = (Y_1000-mean_Y_1000)'*(Y_1000-mean_Y_1000)/length(Y_1000);
var_Y_10000 = (Y_10000-mean_Y_10000)'*(Y_10000-mean_Y_10000)/length(Y_10000);
var_Y_100000 = (Y_100000-mean_Y_100000)'*(Y_100000-mean_Y_100000)/length(Y_100000);

std_dev_Y_100 = sqrt(var_Y_100);
std_dev_Y_1000 = sqrt(var_Y_1000);
std_dev_Y_10000 = sqrt(var_Y_10000);
std_dev_Y_100000 = sqrt(var_Y_100000);

fprintf('Variance and standard deviation for Y (100 trials)\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_Y_100, std_dev_Y_100);
fprintf('\n');
fprintf('Variance and standard deviation for Y (1,000 trials)\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_Y_1000, std_dev_Y_1000);
fprintf('\n');
fprintf('Variance and standard deviation for Y (10,000 trials)\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_Y_10000, std_dev_Y_10000);
fprintf('\n');
fprintf('Variance and standard deviation for Y (100,000 trials)\n');
fprintf('Var: %.2f | Std Dev: %.2f\n', var_Y_100000, std_dev_Y_100000);
fprintf('\n');

% ***********************************************
% 3e

fprintf('\n');
fprintf('Problem 3e *****************************\n');
fprintf('\n');

cdf_100 = [trials_100_probs(1),
           trials_100_probs(1) + trials_100_probs(2),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6) + trials_100_probs(7),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6) + trials_100_probs(7) + trials_100_probs(8),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6) + trials_100_probs(7) + trials_100_probs(8) + trials_100_probs(9),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6) + trials_100_probs(7) + trials_100_probs(8) + trials_100_probs(9) + trials_100_probs(10),
           trials_100_probs(1) + trials_100_probs(2) + trials_100_probs(3) + trials_100_probs(4) + trials_100_probs(5) + trials_100_probs(6) + trials_100_probs(7) + trials_100_probs(8) + trials_100_probs(9) + trials_100_probs(10) + trials_1000_probs(11)];
       
cdf_1000 = [trials_1000_probs(1),
            trials_1000_probs(1) + trials_1000_probs(2),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6) + trials_1000_probs(7),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6) + trials_1000_probs(7) + trials_1000_probs(8),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6) + trials_1000_probs(7) + trials_1000_probs(8) + trials_1000_probs(9),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6) + trials_1000_probs(7) + trials_1000_probs(8) + trials_1000_probs(9) + trials_1000_probs(10),
            trials_1000_probs(1) + trials_1000_probs(2) + trials_1000_probs(3) + trials_1000_probs(4) + trials_1000_probs(5) + trials_1000_probs(6) + trials_1000_probs(7) + trials_1000_probs(8) + trials_1000_probs(9) + trials_1000_probs(10) + trials_1000_probs(11)];
       
cdf_10000 = [trials_10000_probs(1),
             trials_10000_probs(1) + trials_10000_probs(2),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6) + trials_10000_probs(7),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6) + trials_10000_probs(7) + trials_10000_probs(8),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6) + trials_10000_probs(7) + trials_10000_probs(8) + trials_10000_probs(9),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6) + trials_10000_probs(7) + trials_10000_probs(8) + trials_10000_probs(9) + trials_10000_probs(10),
             trials_10000_probs(1) + trials_10000_probs(2) + trials_10000_probs(3) + trials_10000_probs(4) + trials_10000_probs(5) + trials_10000_probs(6) + trials_10000_probs(7) + trials_10000_probs(8) + trials_10000_probs(9) + trials_10000_probs(10) + trials_10000_probs(11)];       

cdf_100000 = [trials_100000_probs(1),
              trials_100000_probs(1) + trials_100000_probs(2),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6) + trials_100000_probs(7),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6) + trials_100000_probs(7) + trials_100000_probs(8),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6) + trials_100000_probs(7) + trials_100000_probs(8) + trials_100000_probs(9),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6) + trials_100000_probs(7) + trials_100000_probs(8) + trials_100000_probs(9) + trials_100000_probs(10),
              trials_100000_probs(1) + trials_100000_probs(2) + trials_100000_probs(3) + trials_100000_probs(4) + trials_100000_probs(5) + trials_100000_probs(6) + trials_100000_probs(7) + trials_100000_probs(8) + trials_100000_probs(9) + trials_100000_probs(10) + trials_100000_probs(11)];

figure5 = figure;
YMatrix1 = [cdf_100,cdf_1000,cdf_10000,cdf_100000];
axes5 = axes('Parent',figure5);
hold(axes5,'on');
plot5 = plot(YMatrix1,'Parent',axes5);
set(plot5(1),'DisplayName','100 trials');
set(plot5(2),'DisplayName','1,000 trials');
set(plot5(3),'DisplayName','10,000 trials');
set(plot5(4),'DisplayName','100,000 trials');
ylabel('Probability');
xlabel('Y value (X1 + X2)');
title('Problem 3e - CDF');
set(axes5,'XTick',[1 2 3 4 5 6 7 8 9 10 11],'XTickLabel',...
    {'2','3','4','5','6','7','8','9','10','11','12'});
legend5 = legend(axes5,'show');
set(legend5,'Location','southeast');

fprintf('see figure 5\n');

%% Problem 4
% ***********************************************
% 4a

fprintf('\n');
fprintf('Problem 4a *****************************\n');
fprintf('\n');

x = .01:.01:100;

lambda1 = 2;
lambda2 = 4;
lambda3 = 10;
lambda4 = 50;

pois1 = poisspdf(x,lambda1)';
pois2 = poisspdf(x,lambda2)';
pois3 = poisspdf(x,lambda3)';
pois4 = poisspdf(x,lambda4)';

YMatrix2 = [pois1,pois2,pois3,pois4];

figure6 = figure;
YMatrix1 = [cdf_100,cdf_1000,cdf_10000,cdf_100000];
axes6 = axes('Parent',figure6);
hold(axes6,'on');
plot6 = plot(YMatrix2,'Parent',axes6);
set(plot6(1),'DisplayName','\lambda=2');
set(plot6(2),'DisplayName','\lambda=4');
set(plot6(3),'DisplayName','\lambda=10');
set(plot6(4),'DisplayName','\lambda=50');
xlabel('k');
ylabel('Probability');
title('Problem 4a');
legend6 = legend(axes6,'show');
set(legend6,'Location','northeast');

fprintf('see figure 6\n');

% ***********************************************
% 4b

fprintf('\n');
fprintf('Problem 4b *****************************\n');
fprintf('\n');

fprintf('λ is both mean and variance for Poisson distros\n');

% ***********************************************
% 4c

fprintf('\n');
fprintf('Problem 4c *****************************\n');
fprintf('\n');

fprintf('Skewness decreases when lambda increases\n');
fprintf('Skewness is defined as γ=1/λ\n');

gamma1 = 1/lambda1;
gamma2 = 1/lambda2;
gamma3 = 1/lambda3;
gamma4 = 1/lambda4;

fprintf('λ=%i,  γ=%.2f\n',lambda1, gamma1);
fprintf('λ=%i,  γ=%.2f\n',lambda2, gamma2);
fprintf('λ=%i, γ=%.2f\n',lambda3, gamma3);
fprintf('λ=%i, γ=%.2f\n',lambda4, gamma4);

% ***********************************************
% 4d

fprintf('\n');
fprintf('Problem 4d *****************************\n');
fprintf('\n');

gaus1 = normpdf(x,lambda1,lambda1)';
gaus2 = normpdf(x,lambda2,lambda2)';
gaus3 = normpdf(x,lambda3,lambda3)';
gaus4 = normpdf(x,lambda4,lambda4)';

YMatrix3 = [pois1,gaus1];
YMatrix4 = [pois2,gaus2];
YMatrix5 = [pois3,gaus3];
YMatrix6 = [pois4,gaus4];
 
figure7 = figure;

subplot1 = subplot(2,2,1,'Parent',figure7);
hold(subplot1,'on');
plot1 = plot(YMatrix3,'Parent',subplot1);
set(plot1(1),'DisplayName','Poisson');
set(plot1(2),'DisplayName','Gaussian');
ylabel({'Probability'});
ylim([0 1]);
xlabel({'k'});
title({'λ=2'});
box(subplot1,'on');
legend(subplot1,'show');


subplot2 = subplot(2,2,2,'Parent',figure7);
hold(subplot2,'on');
plot2 = plot(YMatrix4,'Parent',subplot2);
set(plot2(1),'DisplayName','Poisson');
set(plot2(2),'DisplayName','Gaussian');
ylabel({'Probability'});
ylim([0 1]);
xlabel({'k'});
title({'λ=4'});
box(subplot2,'on');
legend(subplot2,'show');

subplot3 = subplot(2,2,3,'Parent',figure7);
hold(subplot3,'on');
plot3 = plot(YMatrix5,'Parent',subplot3);
set(plot3(1),'DisplayName','Poisson');
set(plot3(2),'DisplayName','Gaussian');
ylabel({'Probability'});
ylim([0 1]);
xlabel({'k'});
title({'λ=10'});
box(subplot3,'on');
legend(subplot3,'show');

subplot4 = subplot(2,2,4,'Parent',figure7);
hold(subplot4,'on');
plot4 = plot(YMatrix6,'Parent',subplot4);
set(plot4(1),'DisplayName','Poisson');
set(plot4(2),'DisplayName','Gaussian');
ylabel({'Probability'});
ylim([0 1]);
xlabel({'k'});
title({'λ=50'});
box(subplot4,'on');
legend(subplot4,'show');
 
fprintf('see figure 7\n');

% ***********************************************
% 4e

fprintf('\n');
fprintf('Problem 4e *****************************\n');
fprintf('\n');

fprintf('Lower lambdas seem to produce better Gaussian approximations\n');

% ***********************************************
% 4f

fprintf('\n');
fprintf('Problem 4f *****************************\n');
fprintf('\n');

figure8 = figure;

subplot1 = subplot(2,2,1,'Parent',figure8);
qqplot(YMatrix3);
title({'λ=2'});
ylabel({'Probability'});
ylim([-0.05 0.3]);

subplot1 = subplot(2,2,2,'Parent',figure8);
qqplot(YMatrix4);
title({'λ=4'});
ylabel({'Probability'});
ylim([-0.05 0.3]);

subplot1 = subplot(2,2,3,'Parent',figure8);
qqplot(YMatrix5);
title({'λ=10'});
ylabel({'Probability'});
ylim([-0.05 0.3]);

subplot1 = subplot(2,2,4,'Parent',figure8);
qqplot(YMatrix6);
title({'λ=50'});
ylabel({'Probability'});
ylim([-0.05 0.3]);

fprintf('see figure 8 for Q-Q plots\n');
fprintf('\n');
fprintf('Higher lambdas are more Gaussian than lower lambda values,\n');
fprintf('demonstrated by the colinearity of the plots\n');

% ***********************************************
% 4g

fprintf('\n');
fprintf('Problem 4g *****************************\n');
fprintf('\n');

fprintf('Poisson distributions converge to a Gaussian distribution as lambda increases\n');