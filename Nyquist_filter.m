function [output] = Nyquist_filter(fs, rate, sample_length, beta)
    
    % filtering
    fc = fs/2;
    T = 1/fs;

    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween

    f_axis = 0:fs/(sample_length*rate):fs - fs/(rate*sample_length);
    f_filter_response = zeros(1, length(f_axis));
    for i = 1:length(f_filter_response)
        if abs(f_axis(i)) <= (1-beta)/(2*T)
            f_filter_response(i) = 1;
        elseif abs(f_axis(i)) > (1-beta)/(2*T) && abs(f_axis(i)) <= (1+ beta)/(2*T)

             f_filter_response(i) = (1*0.5*(1+cos(pi*T/beta * (abs(f_axis(i)) - (1-beta)/(2*T)))));
        else
             f_filter_response(i) = 0;
        end   
    end
    neg_f = fliplr(f_filter_response);
    output = [neg_f f_filter_response];

    figure
    plot([ fliplr(-f_axis)  f_axis], output)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title('Nyquist filter in frequency domain')

    output = f_filter_response

%     t = -T:T/rate:T - T/rate;
%     figure
%     filter_response = ifft(fftshift(output));
%     plot(t , filter_response)
   
end

