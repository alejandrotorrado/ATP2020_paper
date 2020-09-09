function [output_array] = reducezeros(input_array)
%Input array 'input array' to replace all 'NaN' strings with the values in the previous
%rows to produce 'output array'
% keyboard;
a = input_array;
for i = 1:size(a,1)
    for j = 3:size(a,2)
        if a(i,j) == 0
            if i==1
                a(i,j) = a(i+1,j);
            else
                a(i,j) = a(i-1,j);
            end
        end
    end
    
end
output_array = a;
end

