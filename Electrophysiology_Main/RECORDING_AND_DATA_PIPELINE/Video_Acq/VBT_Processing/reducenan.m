function [output_array] = reducenan(input_array)
%Input array 'input array' to replace all 'NaN' strings with the values in the previous
%rows to produce 'output array'
a = input_array;
for i = 1:size(a,1)
    for j = 1:size(a,2)
        if isnan(a(i,j))
            a(i,j) = a(i-1,j);
        end
    end
    
end
output_array = a;
end

