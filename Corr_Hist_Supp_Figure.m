load('./rawData/GLCptspp/case1.mat');
load('Enzymes.mat');

% Create variables name tags and keep 50k samples.
variables = Enzymes;
case1 = samples(:,1:50000)';

tol=10^-9;
% Preprocessing to remove zero variance variables that are 0 in this case.
% We do this as we need as semipositive matrix to run the exact normal
% method. 
index=std(case1)>tol;
variables=variables(index);
data=case1(:,index);

Smat=cov(data); % covariance matrix
Gamma=corrcoef(data); % correlation amtrix

imagesc(Gamma);
set(gca, 'YTick', 1:length(Gamma));
set(gca, 'XTick', 1:length(Gamma));
labelNames = variables;
set(gca,'XTickLabel',labelNames);   % gca gets the current axis
set(gca,'YTickLabel',labelNames);   % gca gets the current axis
set(gca,'XTickLabelRotation',90)
set(gca,'fontsize',2)
red = [1,0,0];
green = [0,1,0];

R = linspace(red(1),green(1),256);
G = linspace(red(2),green(2),256);
B = linspace(red(3),green(3),256);

map = [R', G', B'];
colormap(map)
colorbar

% Plot for PFK distribution as example
idSelected=find(ismember(variables,{'PFK'}));
for i=idSelected%1:50
figure
hist(data(:,i)',100)
title(['C^{GLCptspp}_{',variables{i},'}']);
ylabel('Frequency')
set(gca,'fontsize',22)
set(gca,'fontweight','bold')
grid on
end

% Plot for PFK distribution as example
idSelected=find(ismember(variables,{'PFK','RPE','TPI'}));
for i=idSelected'%1:50
figure
[density, value] = ksdensity(data(:,i),'Kernel','epanechnikov')
plot(value,density,'LineWidth',1.5);
title(['C^{GLCptspp}_{',variables{i},'}']);
ylabel('Density')
set(gca,'fontsize',22)
set(gca,'fontweight','bold')
grid on

end


