function [time, p] = find_lat()
A = 10;
fs = 100;
t = 5;
timeaxis = 0:1/fs:((fs*t)-1)/fs;
% make the data
mu = 5;
sigma = mu*0.20;
y1 = A + sigma*randn(size(timeaxis)) + mu;
sigma = mu*0.36;
y2 = A + sigma*randn(size(timeaxis)) + mu;
y = [y1 y2];

h = zeros(1000,1);
p = zeros(1000,1);

% iterate through data
for i = 101:length(y)
    a = y(i-100:i-51);
    b = y(i-50:i-1);
    [h(i),p(i)] = vartest2(a,b);    % Two-sample F-test
end
% semilogy(p);

% "find" finds any value in array p less than the filter threshold of
% 0.005. It only looks at valid time points (t<1s). Because of that,
% returned values (array indicies) are offset by 100 (sample rate and
% all that).
lat = (find(p(101:1000,:)<0.005)+100);

% latency time *should* be the first index of lat. I could run a check here
% to see if it isn't to test the efficacy of the threshold. Instead I
% just return the first value after 5 seconds
time = lat(find(lat>501,1));