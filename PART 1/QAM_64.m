function [output] = QAM_64(amplitude, input)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% amplitude is max amplitude of the signal
amplitude = abs(amplitude); % little check
if (not(mod(length(input),6) == 0))
    input = [input zeros(1, 6 - mod(length(input), 6))];
end
output = zeros(1, length(input)/6); % add zeros at end of transmission if length of input is uneven

% try to create it using QAM_16

for i = 1:length(input)/6
    if (input(6*i-5) == 0 && input(6*i-2) == 0)
        output(i) = QAM_16(amplitude, [input(6*i - 4: 6*i -3) input(6*i - 1: 6*i)])*3/7 + amplitude*(1/2)*(j-1) + amplitude*(1/14)*(+j-1); %amplitude*(1/2)*(j-1) + amplitude*(1/14)*(+j-1); is offset from QAM_16

    elseif(input(6*i-5) == 0 && input(6*i - 2) == 1)
        temp = QAM_16(amplitude, [input(6*i - 4: 6*i -3) input(6*i - 1: 6*i)])*3/7 + amplitude*(1/2)*(-1-j) + amplitude*(1/14)*(+j-1);
        if(round(imag(temp), 10) == round(-amplitude*6/7, 10))
            output(i) = real(temp) - j*amplitude*1/7;
        else
            output(i) = temp - j*amplitude*3/7;
        end

    elseif(input(6*i-5) == 1 && input(6*i - 2) == 0)
        temp = QAM_16(amplitude, [input(6*i - 4: 6*i -3) input(6*i - 1: 6*i)])*3/7 + amplitude*(1/2)*(1+j) + amplitude*(1/14)*(+j-1);
        if(round(real(temp), 10) == round(amplitude*6/7, 10)) %I need to use round because it is otherwise false
            output(i) = amplitude*1/7 + j*imag(temp);
        else
            output(i) = temp + amplitude*3/7;
        end

   elseif(input(6*i-5) == 1 && input(6*i - 2) == 1)

        [input(6*i - 4: 6*i -3) input(6*i - 1: 6*i)] %debug
        temp = QAM_16(amplitude, [input(6*i - 4: 6*i -3) input(6*i - 1: 6*i)])*3/7 + amplitude*(1/2)*(1-j) + amplitude*(1/14)*(+j-1);
        if(round(real(temp), 10) == round(amplitude*6/7, 10))
            output(i) = amplitude*1/7;
        else
            output(i) = real(temp) + amplitude*3/7;
        end

        if(round(imag(temp), 10) == round(-amplitude*6/7, 10))
            output(i) = output(i) - j*amplitude*1/7;
        else
            output(i) = output(i) + j*imag(temp) - j*amplitude*3/7;
        end
    else
        output = [];
        break
    end
end