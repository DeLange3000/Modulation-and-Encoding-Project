function [output] = QPSK(amplitude, input)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
amplitude = abs(amplitude); % little check
input = [input zeros(1, mod(length(input), 2))] % add zeros at end of transmission if length of input is uneven
output = zeros(1, length(input)/2);

for i = 1:length(input)/2
    if (input(i*2-1) == 1 && input(2*i) == 1)
        output(i) = -amplitude;
    elseif (input(2*i-1) == 0 && input(2*i) == 0)    
        output(i) = amplitude;
    elseif (input(2*i-1) == 1 && input(2*i) == 0)    
        output(i) = - j*amplitude;
    elseif (input(2*i-1) == 0 && input(2*i) == 1)    
        output(i) = j*amplitude;
    else
        output = [];
        break
    end
end

end