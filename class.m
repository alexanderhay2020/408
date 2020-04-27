% sampsizepwr
% 5 arguments, enter 4 and leave 1 blank

% paired t-test with mean of 100 and std of 5
% you want a power of 8.9 for one-tailed test
% to distinguish mean of 102 from 100
n = sampsizepwr('t', [100 5], 102, .98, [], 'tail', 'right')
power = sampsizepwr('t', [100 5], 102, [], 28,  'tail', 'right')
p1_mean = sampsizepwr('t', [100 5], [], .9, 20,  'tail', 'right')

