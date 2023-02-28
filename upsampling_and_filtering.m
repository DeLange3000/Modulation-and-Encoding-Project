function [output] = upsampling_and_filtering(input, rate)

    % input: complex symbols from modulation
    % oversampling: 1M symbols/s -> in Matlab: use 10 samples per symbol
    oversampled = zeros(1, rate*length(input));
    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;
    for i = 1:length(input)
        oversampled(rate*(i-1)+1:rate*i) =  input(i);
    end

    t = 1/(rate*fs):1/(rate*fs):(length(oversampled))/(rate*fs);
    figure
    plot(t, oversampled)
    title('time domain signal (not filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')


    f_oversampled = (fft(oversampled, length(oversampled)));
    f_axis = 0:rate*fs/length(f_oversampled):rate*fs-rate*fs/length(f_oversampled);
    figure
    plot(f_axis, abs(f_oversampled))
    title('frequency domain signal (not filtered)')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')

    

    h = rcosdesign(0.3,1,100,"sqrt");
    figure
    freqz(h)
    filtered = upfirdn(input, h, rate);
    f_filtered = fft(filtered);


    figure
    hold on
    plot(0:rate*fs/length(f_oversampled):rate*fs, abs(fftshift(f_filtered)))
    title('frequency domain signal (filtered)')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    
    figure
    plot(0:T/rate:length(f_filtered)*T/rate - T/rate, abs(ifft(f_filtered)))
    title('time domain signal (filtered)')
    xlabel('Time (s)')
    ylabel('Amplitude')

    output = ifft(f_filtered);

end