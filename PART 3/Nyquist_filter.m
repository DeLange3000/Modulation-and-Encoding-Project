function [output] = Nyquist_filter(fs, rate, sample_length, beta)

    %% filter parameters
    fc = fs/2;
    T = 1/fs;

	filter_res = 101*rate;

    %% filter calculation

    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween

    f_axis = linspace(-fs*rate,fs*rate, filter_res);
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

      filter_response = ifft(sqrt(ifftshift(f_filter_response)));
      filter_response = fftshift(filter_response); % normalize
      filter_response = filter_response/norm(filter_response);
	  %filter_response = downsample(filter_response,	100)
     figure
     plot((filter_response))
     xlabel('Time (s)')
     ylabel('Amplitude')
     title('Nyquist filter in time domain')
     %figure
     %impz(filter_response)

    %% return positive frequencies of filter

    %output = f_filter_response(length(f_filter_response)/2 + 1:length(f_filter_response))';
    output = filter_response';

end

