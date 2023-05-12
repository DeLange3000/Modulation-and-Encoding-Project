function cfo_est = cfo_est(y, a, argmaxn, K, T)
sum = 0;
for k = 1:K
	sum = sum + angle(diffCorr(y, a, k, argmaxn))/k;
end
cfo_est = -sum/(2*pi*K*T);
end