function output = diffCorr(y, a, k, n)
	N = length(a);
	sum  = 0;
	for l = k+1:N
		y1 = y(n+l);
		y2 = y(n+l-k);
		sum = sum +	conj(y1)*a(l) * conj(conj(y2)*a(l-k));
	end
	output = sum/(N-k);
end