function [output] = filtering_and_downsampling(input,rate, f_filter)

    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;

    %% filter received signal

    f_input = (fft(input));
    f_axis = 0:rate*fs/length(input):rate*fs-rate*fs/length(f_input);
    t = 1/(rate*fs):1/(rate*fs):(length(input))/(rate*fs);
   
    f_filtered = f_input.*f_filter;

    %% plot filtered signal

%     figure
%     hold on
%     plot(t , real(ifft(f_filtered)))
%     %plot(t , imag(ifft(f_filtered)))
%     title('time domain signal (filtered at receiver)')
%     xlabel('Time (s)')
%     ylabel('Amplitude')

    %% downsample

    filtered = (ifft(f_filtered));

%     output = [];
%     for i =1:rate:length(filtered)
%         output = [output; filtered(i)];
%     end

    output = zeros(length(filtered)/rate, 1);
    a = 1;
    for i =1:rate:length(filtered)
        output(a) = filtered(i);
        a = a + 1;
    end
