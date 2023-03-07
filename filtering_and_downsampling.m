function [output] = filtering_and_downsampling(input,rate, f_filter)

    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;

    f_input = (fft(input, length(input)));
    f_axis = 0:rate*fs/length(input):rate*fs-rate*fs/length(f_input);
    t = 1/(rate*fs):1/(rate*fs):(length(input))/(rate*fs);
   


    f_filtered = f_input.*f_filter;


    %plotting
    figure
    plot(t , real(ifft(f_filtered)))
    title('time domain signal (filtered at receiver)')
    xlabel('Time (s)')
    ylabel('Amplitude')

filtered = ifft(f_filtered)
output = zeros(1, length(filtered)/100);
output(1) = filtered(1);
for i =100:100:length(filtered)
    output(i/100) = filtered(i);

end
