function [output] = Nyquist_filter(fs, rate, sample_length, beta, taps)

    %% filter parameters
    fc = fs/2;
    T = 1/(fs);

    %% filter calculation

    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween
	
	lower = -(taps-1)/2;
    f_axis = (lower:-lower)*(fs*rate)/(taps);
    f_filter_response = zeros(1, length(f_axis));
    for i = 1:length(f_filter_response)
        if abs(f_axis(i)) <= (1-beta)/(2*T)
            f_filter_response(i) = T;
        elseif abs(f_axis(i)) > (1-beta)/(2*T) && abs(f_axis(i)) <= (1+ beta)/(2*T)

             f_filter_response(i) = (T*0.5*(1+cos(pi*T/beta * (abs(f_axis(i)) - (1-beta)/(2*T))))); %(1/(T*fc)*0.5*(1+cos(pi*T/beta * (abs(f_axis(i)) - (1-beta)/(2*T)))));
        else
             f_filter_response(i) = 0;
        end   
    end
    
    % plot filter in frequency domain  
    
    figure
    plot(f_axis, f_filter_response)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title('Nyquist filter in frequency domain')

    %% plot impulse response of filter

	t_axis = (lower:-lower)*T/rate;

     f_r = ifft((ifftshift(f_filter_response)));
     filter_response = ifft(sqrt(ifftshift(f_filter_response)));
     filter_response = fftshift(filter_response/sqrt(max(f_r))); % normalize
     figure
     plot(t_axis, (filter_response))
     xlabel('Time (s)')
     ylabel('Amplitude')
     title('Nyquist filter in time domain')
     figure
     impz(filter_response)

    %% return positive frequencies of filter

    %output = f_filter_response(length(f_filter_response)/2 + 1:length(f_filter_response))';
    output = filter_response';

end

