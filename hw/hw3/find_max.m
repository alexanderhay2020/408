function [max_val, i_index, j_index] = find_max(array)
[x y] = size(array);
% x = temp(1);
% y = temp(2);
max_val = -inf;
for i = 1:x
    for j = 1:y
        if array(i,j) > max_val;
            max_val = array(i,j);
            i_index = i;
            j_index = j;
        end
    end
end
