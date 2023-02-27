function output = Nyquist(input, rate)
    % input: complex symbols from modulation
    % oversampling: 1M symbols/s -> in Matlab: use 10 samples per symbol
    oversampled = zeros(1, rate*length(input));
    fc = 1e6;
    for i = 1:length(input)
        oversampled(rate*(i-1)+1:rate*i) = real(input(i)*exp(j*2*pi*fc*(0:1/(rate*fc):(rate-1)/(rate*fc))));
    end
    

    t = 0:1/(rate*fc):(length(oversampled)-1)/(rate*fc);
    plot(t, oversampled)


    f_oversampled = (fft(oversampled, length(oversampled)))
    f_axis = 0:rate*fc/length(f_oversampled):rate*fc-rate*fc/length(f_oversampled);
    figure
    plot(abs(f_axis), abs(f_oversampled))


    % filtering
    beta = 0.3;
    
    

end