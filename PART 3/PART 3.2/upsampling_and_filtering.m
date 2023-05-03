function [output] = upsampling_and_filtering(input, rate, filter)

    %% upsampling

    % input: complex symbols from modulation
    oversampled = zeros(rate*length(input), 1);
    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;
    for i = 1:length(input)
        oversampled(rate*(i-1)+1) = (input(i));
        oversampled(rate*(i-1)+2:rate*(i)) = 0; %input(i);
    end

    oversampled = [zeros(rate - 1, 1); oversampled];

    %% plot unfiltered signal in time and frequency domain

    t = 1/(rate*fs):1/(rate*fs):(length(oversampled))/(rate*fs);
    figure
    hold on
    plot(t, real(oversampled))
    %plot(t, imag(oversampled))
    title('time domain signal (not filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')


    %% filter signal

    filtered = conv(oversampled, filter, 'same');


    %% plot filtered signal
    
    figure
    hold on
    plot(t, real(filtered))
    %plot(t, imag(fftshift(ifft(f_filtered))))
    %plot(0:T/rate:length(f_filtered)*T/rate - T/rate, abs(ifft(f_filtered)))
    title('time domain signal (filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')

%     figure
%     hold on
%     plot( -fs:fs*2/(length(input)*rate):fs - fs/(rate*length(input)), real(fftshift(fft(filtered))))
%     title('frequency domain signal (filtered)')
%     xlabel('Frequency (f)')
%     ylabel('Amplitude')


    %% return filtered signal

    output = filtered;

end