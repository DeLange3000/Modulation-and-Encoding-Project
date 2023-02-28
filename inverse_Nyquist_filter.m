function [output] = inverse_Nyquist_filter(input, rate)

    fc = 1e6;
    f_input = (fft(input, length(input)));
    f_axis = 0:rate*fc/length(input):rate*fc-rate*fc/length(f_input);
    beta = 0.3;

    filter = ones(1, length(input)/2);
    for i = 1:length(f_axis)/2
        if f_axis(i) > (1+beta)*fc/2
            % set to zero
            filter(i) = 0;
        elseif f_axis(i) > (1-beta)*fc/2

            filter(i) = sqrt(0.5*(1+cos(pi/(beta*fc) * (f_axis(i) - (1-beta)*fc/2))));
        end
    end

    % some tweaking for symmetry
    filter = [1 filter(1:end-1) fliplr(filter)];

    %filter
    f_filtered = f_input .* filter;

    %plotting
    t = 1/(rate*fc):1/(rate*fc):(length(input))/(rate*fc);
    figure
    plot(t, input)
    title('time domain signal (filtered at receiver)')
    xlabel('Time (s)')
    ylabel('Amplitude')

    figure
    plot(f_axis, abs(fftshift(f_filtered)))
    title('frequency domain signal (filtered at receiver)')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')

    filtered = ifft(f_filtered);

    %downsampling
    
    output = zeros(1,length(filtered)/rate)
    output(1) = abs(filtered(1)) + j *(phase(filtered(1)));
    for i = 100:100:length(filtered)
            output(i/100) = abs(filtered(i)) + j *(phase(filtered(i))); % gets real part of signal
    end

    output = output';


end