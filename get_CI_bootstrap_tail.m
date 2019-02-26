function CIs = get_CI_bootstrap_tail(data,nsimu,alpha)
% data = [samples x variables] = [n x p];
% nsimu = number of bootsrap resamplings. If not provided use as many as we have samples.
% alpha = confidence level. 5% is defaults in stats. 

if nargin < 3
    alpha=0.05;
end

[n,p]=size(data);

if nargin < 2
   nsimu = n; % num simulations  
end

data_mean=mean(data);

Smat=cov(data); % covariance matrix

t_stat = zeros(nsimu,p);
for i = 1:nsimu
    index=randsample(n,n,true);
    x_star=data(index,:);
    x_star_mean=mean(x_star);
    x_star_std=std(x_star)/sqrt(n);
    t_stat(i,:)=(x_star_mean-data_mean)./x_star_std; % no abs here
end

h_boot=[tiedrank(t_stat)-1]/[nsimu+1];
h_max=max(abs(h_boot-0.5)'); %take abs here instead

q_alpha=quantile(h_max,1-alpha);
beta = 1-2*q_alpha;

u_vec=quantile(t_stat,1-beta/2,1);
l_vec=quantile(t_stat,beta/2,1); % Top values become negative....

Low_tail = data_mean - sqrt(diag(Smat)/n)'.*u_vec;
Upp_tail = data_mean - sqrt(diag(Smat)/n)'.*l_vec;

CIs = [Low_tail; Upp_tail];

end