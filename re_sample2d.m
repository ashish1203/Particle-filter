function p = re_sample2d(xp,pdf) % resampling
	cdf = cumsum(pdf); % cumulative sum
	diff = cdf'*ones(1,length(pdf)) - ones(length(pdf),1)*rand(1,length(pdf));
	diff = (diff <= 0) * 2 + diff;
	[~, idx] = min(diff);
	p = xp(:,idx);
end