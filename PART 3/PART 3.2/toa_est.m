function [argmaxn, cfo_est] = toa_est(a, y, K, T)
sums = zeros(length(y)-length(a), 1);
for n = 1:length(y)-length(a)
sum = 0;
for k = 1:K
	sum = sum + abs(diffCorr(y, a, k, n));
end
sums(n) = sum;
end

stem(sums)
[~, argmaxn] = max(sums);

sum = 0;
for k = 1:K
	sum = sum + angle(diffCorr(y, a, k, argmaxn))/k;
end
cfo_est = -sum/(2*pi*K*T);
end