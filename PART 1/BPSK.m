function [output] = BPSK(amplitude, input)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = zeros(1, length(input));

for i = 1:length(input)
    if (input(i) == 1)
        output(i) = abs(amplitude);
    elseif (input(i) == 0)
        output(i) = -abs(amplitude);
    else
        output = [];
        break
    end
end

end