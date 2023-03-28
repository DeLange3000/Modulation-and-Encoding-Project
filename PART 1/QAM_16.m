function [output] = QAM_16(amplitude,input)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% amplitude is max amplitude of the signal
amplitude = abs(amplitude); % little check
if (not(mod(length(input),4) == 0))
    input = [input zeros(1, 4 - mod(length(input), 4))];
end
output = zeros(1, length(input)/4); % add zeros at end of transmission if length of input is uneven


for i = 1:length(input)/4
    if(input(4*i-3:4*i-2) == [0 0])
        output(i) = -amplitude;
    elseif(input(4*i-3:4*i-2) == [0 1])
        output(i) = -amplitude/3;
    elseif(input(4*i-3:4*i-2) == [1 1])
        output(i) = amplitude/3;
    elseif(input(4*i-3:4*i-2) == [1 0])
        output(i) = amplitude;
    else
        output = [];
        break
    end

    if(input(4*i-1:4*i) == [0 0])
        output(i) = output(i) + j*amplitude;
    elseif(input(4*i-1:4*i) == [0 1])
        output(i) = output(i) + j*amplitude/3;
    elseif(input(4*i-1:4*i) == [1 1])
        output(i) = output(i) - j*amplitude/3;
    elseif(input(4*i-1:4*i) == [1 0])
        output(i) = output(i) - j*amplitude;
    else
        output = [];
        break
    end
end
end