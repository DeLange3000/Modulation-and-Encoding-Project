function [output] = Nyquist_filter(input, rate)
    % input: complex symbols from modulation
    % oversampling: 1M symbols/s -> in Matlab: use 10 samples per symbol
    oversampled = zeros(1, rate*length(input));
    fc = 1e6;
    for i = 1:length(input)
        oversampled(rate*(i-1)+1:rate*i) =  real(input(i)*exp(j*2*pi*fc*(0:1/(rate*fc):(rate-1)/(rate*fc))));
    end

    t = 1/(rate*fc):1/(rate*fc):(length(oversampled))/(rate*fc);
    figure
    plot(t, oversampled)
    title('time domain signal (not filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')


    f_oversampled = (fft(oversampled, length(oversampled)));
    f_axis = 0:rate*fc/length(f_oversampled):rate*fc-rate*fc/length(f_oversampled);
    figure
    plot(f_axis, abs(f_oversampled))
    title('frequency domain signal (not filtered)')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')

    
    % filtering
    beta = 0.3;
    fc = 2e6;

    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween

    % building the filter: first make one half, then mirror it
    filter = ones(1, length(f_oversampled)/2);
    for i = 1:length(f_axis)/2
        if f_axis(i) > (1+beta)*fc/2
            % set to zero
            filter(i) = 0;
        elseif f_axis(i) > (1-beta)*fc/2

            filter(i) = (0.5*(1+cos(pi/(beta*fc) * (f_axis(i) - (1-beta)*fc/2))));
        end
    end

    % some tweaking for symmetry
    filter = [1 filter(1:end-1) fliplr(filter)];
    figure
    subplot(221)
    plot(f_axis, filter)
    title("frequency resp")
    subplot(222)
    plot(t, ifftshift(ifft(filter)))
    title("impulse response")
    subplot(223)
    title("")


    %f_filtered = f_oversampled .* filter;
    %filtered = (ifft(f_filtered, length(oversampled)));

    f
    figure
    plot(fftshift(abs(f_oversampled)))
    hold on
    plot(abs(fftshift(f_filtered)))
    plot(fftshift(filter))
    title('frequency domain signal (filtered)')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    legend('unfiltered signal', 'filtered signal')
    
    figure
    plot(t, filtered)
    title('time domain signal (filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')

    output = filtered;
    

end