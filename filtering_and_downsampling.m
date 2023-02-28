function [output] = filtering_and_downsampling(input,rate)

    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;

    f_input = (fft(input, length(input)));
    f_axis = 0:rate*fs/length(input):rate*fs-rate*fs/length(f_input);
   

    h = rcosdesign(0.3,1,100,"sqrt");
    filtered = upfirdn(input, h, 1, rate);


    %plotting
    figure
    plot(0:T/rate:length(filtered)*T/rate - T/rate , filtered)
    title('time domain signal (filtered at receiver)')
    xlabel('Time (s)')
    ylabel('Amplitude')


output = filtered(2:length(filtered)-1);


end
