function [output, Garner_errors] = filtering_and_downsampling(input,rate, filter, shifts, error_factor_K)

    fc = 1e6;
    fs = 2e6; %symbol rate
    T = 1/fs;
    Tsample = 1/(rate*fs);
    Garner_error = 0;
    Garner_errors = zeros(1,length((rate + shifts - Garner_error):rate:length(input)));

    %% filter received signal

    t = 1/(rate*fs):1/(rate*fs):(length(input))/(rate*fs);
   
    filtered = conv(input, fliplr(filter), 'same');

    %% plot filtered signal
% 
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


    output = zeros((length(filtered)-(rate-1))/rate, 1);
    a = 1;
    Garner_error_index = 0;
    for i = (rate + shifts):rate:length(filtered)
        i = i + Garner_error_index;
        if(i > length(filtered))
            uh = 0;
        end
        %plot(t(round(i)), real(filtered(round(i))), '*')
        if(mod(i, 1) ~= 0)
            %interpolation
            i1 = floor(i);
            i2 = ceil(i);
            output(a) = (filtered(i1)*(i2 - i)+filtered(i2)*(i - i1))/(i2 - i1);
        else
            output(a) = filtered(i);
           
        end

        if(i == rate + shifts)
            Garner_error = 0;
        else
            %get sample halfway
            if(mod(i - rate/2, 1) ~= 0)
                %interpolation
                i1 = floor(i - rate/2);
                i2 = ceil(i - rate/2);
                sample_halfway = (filtered(i1)*(i2 - i + rate/2)+filtered(i2)*(i - rate/2 - i1))/(i2 - i1);
            else
                sample_halfway = filtered(i - rate/2);
            end
            Garner_error = Garner_error - 2*error_factor_K*real(sample_halfway*(conj(output(a))- conj(output(a-1))));
            Garner_errors(a) =  Garner_error;
            Garner_error_index = Garner_error/t(1)*T;

        end
        a = a + 1;

    end
end
