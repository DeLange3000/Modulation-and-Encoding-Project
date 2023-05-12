function [argmaxn, cfo] = toa_est(y, a, K, T)
sums = zeros(length(y)-length(a), 1);
for n = 1:(length(y)-length(a))
sum = 0;
for k = 1:K
	sum = sum + abs(diffCorr(y, a, k, n));
end
sums(n) = sum;
end
figure
stem(1:(length(y)-length(a)), sums)
[~, argmaxn] = max(sums);

cfo = cfo_est(y, a, argmaxn, K, T);
end