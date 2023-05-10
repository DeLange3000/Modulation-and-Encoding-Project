a = 2*(randi(2, 1, 10)-1.5);
y = [2*(randi(2, 1, 100)-1.5) a 2*(randi(2, 1, 500)-1.5)];
K = length(a)-1;
%[n, cfo] = toa_est(a, y, K);
clc

pilot = [0 0 0 0];
pilot_length = length(pilot);
pilot_int = 7;
bitstream = ones(1, 10);
bitstream_length = length(bitstream);



for i = 1:pilot_int:bitstream_length+pilot_length*ceil(bitstream_length/(pilot_int-pilot_length))
	i
	begin = bitstream(1:i-1)
	rest = bitstream(i:end)
	bitstream = [begin pilot rest]
end

bitstream_length+pilot_length*ceil(bitstream_length/(pilot_int-pilot_length))
length(bitstream)