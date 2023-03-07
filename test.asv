clear
close all
clc

fs = 2e6
T = 1/fs

f_signal = zeros(1, 1000);
f_signal(1) = 1;
f = 0:fs/1000:fs-fs/1000;

f_filtered_signal = zeros(1, length(f));
for i = 1:length(f)
    f_filtered_signal(i) = Nyquist_filter(f_signal(i), f(i), 100);
end

plot(real(ifft(f_filtered_signal)))

figure
plot(f, abs(f_filtered_signal))



