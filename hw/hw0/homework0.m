% Alexander Hay
% Homework 0

%% Question 1
% installed

%% Question 2
x=6;        % this stores the variable
x=6         % this prints to the command window
y=2;   

x+y;        % ans=8
x-y;        % ans=4
x*y;        % ans=12
x/y;        % and=3

% variables are x, y, and ans

%% Question 3
a=[1,2,5,7,4];
a(1);       % returns 1, first index in vector
%a(0);      % returns an error, matlab indicies start at 1
a([2,4]);   % returns 2 and 7, the 2nd and 4th indicies of the vector
a+y;        % adds y to each value in vector, returns vector
            % [3,4,7,9,6]
a*y;        % multiplies y by each value in vector, returns vector
            % [2,4,10,14,8]
length(a);  % returns 5, number of elements in vector
min(a);     % returns 1
max(a);     % returns 7

%% Question 4
A = [1,2,3;4,5,6];
A(1);       % returns 1
A(1,1);     % returns 1
A(1,2);     % returns 2
A(2,1);     % returns 4
A*y;        % returns 2x3 matrix
            % [2, 4, 6]
            % [8,10,12]
A'          % returns A transposed