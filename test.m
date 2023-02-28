

filter = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          0.3, ...
  'FilterSpanInSymbols',    1, ...
  'OutputSamplesPerSymbol', 100);

a = filter.coeffs()

filter2 = rcosdesign(0.3, 10, 100, 'sqrt');

output = upfirdn(encoded_signal, filter2, 100);

figure
plot(abs(output))
figure
plot(abs(fft(output)))
