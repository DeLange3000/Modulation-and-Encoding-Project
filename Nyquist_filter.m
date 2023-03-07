function [output] = Nyquist_filter(input, f_sample, rate)
    
    % filtering
    beta = 0.3;
    fc = 1e6;
    fs = 2e6;
    T = 1/fs;

    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween

    % building the filter: first make one half, then mirror it
        if abs(f_sample) > (1+beta)*fs/2
            % set to zero
            output = 0;
        elseif abs(f_sample) > (1-beta)*fs/2

            output = sqrt((T*0.5*(1+cos(pi/(beta*fs) * (abs(f_sample) - (1-beta)*fs/2)))))*input;
        else
            output = sqrt(T);
        end   
end