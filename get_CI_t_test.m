function CIs = get_CI_t_test(data,alpha)
% data = [samples x variables] = [n x p];
% alpha = confidence level. 5% is the default.

if nargin < 2
    alpha=0.05;
end

[n,p]=size(data);

data_mean=mean(data);

CIs = repmat(data_mean,2,1) + [-1; 1]*tinv(1-alpha/2,n-1)...
    *[std(data)/sqrt(n)];

end