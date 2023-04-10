function [output] = filtering_and_downsampling(input,rate, filter)

    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;

    %% filter received signal

    t = 1/(rate*fs):1/(rate*fs):(length(input))/(rate*fs);
   
    filtered = conv(input, fliplr(filter), 'same');

    %% plot filtered signal

%     figure
%     hold on
%     plot(t , real(filtered))
%     %plot(t , imag(ifft(f_filtered)))
%     title('time domain signal (filtered at receiver)')
%     xlabel('Time (s)')
%     ylabel('Amplitude')

    %% downsample

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
