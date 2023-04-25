function filter = rcfilt(beta, fc, f_axis)
    % freq resp of filter:
    % 1 for frequencies up to (1-beta)fc/2
    % 0 for frequencies from (1+beta)fc/2 to infinity
    % sqrt(0.5*(1+cos(pi/(beta*fc) * (f - (1-beta)fc/2))) for frequencies
    % inbetween

    % building the filter: first make one half, then mirror it
    filter = ones(1, length(f_axis)/2);
    for i = 1:length(f_axis)/2
        if f_axis(i) > (1+beta)*fc/2
            % set to zero
            filter(i) = 0;
        elseif f_axis(i) > (1-beta)*fc/2

            filter(i) = (0.5*(1+cos(pi/(beta*fc) * (f_axis(i) - (1-beta)*fc/2))));
        end
    end

    % some tweaking for symmetry
    filter = [ filter(1:end) fliplr(filter)];

    figure
    plot(f_axis, abs(filter))
    
end