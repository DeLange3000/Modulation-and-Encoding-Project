clc;
close all

% random bitsstream generator
bit_stream = [];
for i = 1:30
    bit_stream = [bit_stream round(rand())];
end

%bit_stream = [0 0 1 1 1 1]
%output = BPSK(1, bit_stream);
%output = QPSK(1, bit_stream);
%output = QAM_16(1, bit_stream);
output = QAM_64(1, bit_stream);


if (length(output) == 0)
    fprintf('an error has occured')
else
%     figure
%     hold on
%     stem(real(output));
%     stem(imag(output));
    figure
    plot(real(output), imag(output), '*');
    title('Modulated signal')
    xlabel('Real axis')
    ylabel('Imaginairy axis')
    set(gca, 'XAxisLocation', 'origin')
    set(gca, 'YAxisLocation', 'origin')
end

filtered_signal = Nyquist_filter(output,100)

noisy_signal = Add_noise(filtered_signal)



