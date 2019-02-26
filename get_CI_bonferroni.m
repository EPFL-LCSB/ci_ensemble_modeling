function CIs = get_CI_bonferroni(data,alpha)
% data = [samples x variables] = [n x p];
% alpha = confidence level. 5% is defaults in stats. 

if nargin < 2
    alpha=0.05;
end

[n,p]=size(data);

data_mean=mean(data);

sd_vec = std(data);
alpha_bonf = alpha/p;
Z_val=norminv(1-[alpha_bonf/2]);
% Print num samples required for 0.1 EM 
disp('Samples needed for EM of 0.1 !!')
max((sd_vec*Z_val/0.1).^2) % 0.1 here is the error margin EM. 

CIs = repmat(data_mean,2,1) + [-1; 1]*tinv(1-[alpha/p/2],n-1)...
    *[std(data)/sqrt(n)];

end