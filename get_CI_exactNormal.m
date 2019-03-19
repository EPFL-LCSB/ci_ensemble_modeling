function CIs = get_CI_exactNormal(data,nsimu,alpha)
% data = [samples x variables] = [n x p];
% nsimu = number of simulations using normal distribution with implied
% variance matrix. If not provided use as many as we have samples.
% alpha = confidence level. 5% is defaults in stats. 

% REMOVE ZERO STD VARIABLES AS mvnrnd WILL THROW AN ERROR FOR NOT BEING
% ABLE TO GENERATE POSITIVE SEMI DEFINITE MATRIX

if nargin < 3
    alpha=0.05;
end

[n,p]=size(data);

if nargin < 2
   nsimu = n; % num simulations  
end

data_mean=mean(data);

Smat=cov(data); % covariance matrix
Gamma=corrcoef(data); % correlation matrix
Zsimu = mvnrnd(zeros(p,1),Gamma,nsimu);% simulate Z nsimu times
absmax = max(abs(Zsimu)'); %compute max_j |Z_j|
q_alpha=quantile(absmax,1-alpha); 

Low = data_mean - q_alpha*sqrt(diag(Smat)/n)';
Upp = data_mean + q_alpha*sqrt(diag(Smat)/n)';

CIs = [Low; Upp];

end