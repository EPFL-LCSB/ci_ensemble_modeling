
% Set the random seed
RandStream.setGlobalStream(RandStream('mt19937ar','seed',2011));

% load data and corresponding enyzme name tags
load('./rawData/GLCptspp/case1.mat');
load('Enzymes.mat');

% Create variables name tags and keep all samples. Could use less to make
% run bit fasterto test.
variables = Enzymes; % enzymatic reactions
case1 = samples(:,1:50000)'; % [samples x vars]

% Preprocessing to remove zero variance variables that are 0 in this case.
% We do this as we need as semipositive matrix to run the exact normal
% method.
tol=10^-9;
index=std(case1)>tol;
variables=variables(index);
case1=case1(:,index);

% clear variables not required
clear Enzymes samples

% Initiate varible to collect confidence intervals.
CIagg=[];

% Classical univariate confidence intervals (CI) without any correction for
% simulataneous comparison. 95% confidence level.
% Note: with multiple variables the t-distribution tends towards the normal
% distribution.
CIagg.ttest=get_CI_t_test(case1);

% Bonferroni correction of the t-test CIs for simultaneous testing.
startBonf = tic;
CIagg.bonf=get_CI_bonferroni(case1);
tBonf = toc(startBonf);

% Exact normal method for correcting the CIs. Here we account for the
% interdependencies of the variables witht he correlation matrix. We still
% make the normalty assumption though.
startExact = tic;
CIagg.norm=get_CI_exactNormal(case1);
tExact = toc(startExact);

% Bootstrappping approach for constructing CIs. We do not rely on normality
% assumption anymore as the CIs are built around a "pivot". A pseudo
% "t-statistic" is computed from the provided data to estiamte the CIs via
% resampling from the data with replacement.
startBoot = tic;
CIagg.boot=get_CI_bootstrap_tail(case1,25000); 
tBoot = toc(startBoot);

% data mean
data_mean = mean(case1);

% Plot some of these CIs
plot_CI(data_mean,CIagg.ttest)
plot_CI(data_mean,CIagg.bonf)
plot_CI(data_mean,CIagg.norm)
plot_CI(data_mean,CIagg.boot)

% Select top candidates based on absolute magnitude for plotting and
% comparison.
numTops=10;
[~,m_vec_order] = sort(abs(data_mean),'descend');
top_ID = m_vec_order(1:numTops);

% Compare the three methods
figure
errorbar([1:3:numTops*3],data_mean(top_ID),data_mean(top_ID)-CIagg.bonf(1,top_ID)...
    ,CIagg.bonf(2,top_ID)-data_mean(top_ID),'db','MarkerFaceColor','b','MarkerSize',8,'LineWidth',1.5);
hold on
errorbar([2:3:numTops*3],data_mean(top_ID),data_mean(top_ID)-CIagg.norm(1,top_ID)...
    ,CIagg.norm(2,top_ID)-data_mean(top_ID),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1.5);
errorbar([3:3:numTops*3],data_mean(top_ID),data_mean(top_ID)-CIagg.boot(1,top_ID)...
    ,CIagg.boot(2,top_ID)-data_mean(top_ID),'dk','MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5);
grid on;

set(gca,'XTick',[2:3:numTops*3])
set(gca,'XTickLabel',variables(top_ID))
set(gca,'XTickLabelRotation',45)
plot(xlim,[0 0],'--k')
ylabel('Means and CIs')
set(gca,'fontsize',22)
set(gca,'fontweight','bold')

legend('Bonferroni','Exact normal','Bootstrap')

% Alternative plot - compare the three methods to the univariate no correction
figure
errorbar([1:4:numTops*4],data_mean(top_ID),data_mean(top_ID)-CIagg.ttest(1,top_ID)...
    ,CIagg.ttest(2,top_ID)-data_mean(top_ID),'dm','MarkerFaceColor','m','MarkerSize',8,'LineWidth',1.5);
hold on
errorbar([2:4:numTops*4],data_mean(top_ID),data_mean(top_ID)-CIagg.bonf(1,top_ID)...
    ,CIagg.bonf(2,top_ID)-data_mean(top_ID),'db','MarkerFaceColor','b','MarkerSize',8,'LineWidth',1.5);
errorbar([3:4:numTops*4],data_mean(top_ID),data_mean(top_ID)-CIagg.norm(1,top_ID)...
    ,CIagg.norm(2,top_ID)-data_mean(top_ID),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1.5);
errorbar([4:4:numTops*4],data_mean(top_ID),data_mean(top_ID)-CIagg.boot(1,top_ID)...
    ,CIagg.boot(2,top_ID)-data_mean(top_ID),'dk','MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5);
grid on;

set(gca,'XTick',[2.5:4:numTops*4])
set(gca,'XTickLabel',variables(top_ID))
set(gca,'XTickLabelRotation',45)
plot(xlim,[0 0],'--k')
ylabel('Means and CIs')
set(gca,'fontsize',22)
set(gca,'fontweight','bold')

legend('Univariate','Bonferroni','Exact normal','Bootstrap')

%% Case studiy: applying the three statistical methods
% for constructing CIs when comparing 4 diffferent FDPs.
% We select top 7 enzymes per case and take the union of these top enzymes for study.

% Load the names to fetch variable names
load('Enzymes.mat');

% Create variables name tags and keep 50k samples.
variables = Enzymes;
noSamples=50000;
load('./rawData/GLCptspp/case1.mat');
case1 = samples(:,1:noSamples)'; % [samples x vars]
load('./rawData/GLCptspp/case2.mat');
case2 = samples(:,1:noSamples)'; % [samples x vars]
load('./rawData/GLCptspp/case3.mat');
case3 = samples(:,1:noSamples)'; % [samples x vars]
load('./rawData/GLCptspp/case4.mat');
case4 = samples(:,1:noSamples)'; % [samples x vars]

% clear variables not required
clear Enzymes samples

% Find top variables using the bootstrapping approach. We favour this
% method as the distribution of control coefficients are generally not very
% normal...
noTopVar=7;
tol=10^-9;
for i=1:4
    % evaluate eache case to get significant variables
    eval(['index=std(case',num2str(i),')>tol;']) % consider ones that present variance larger than tolerance...
    eval(['dat = case',num2str(i),'(:,index);'])
    
    m_vec = mean(dat); % get means
    
    [~,m_vec_order]=sort(abs(m_vec),'descend'); % order absolute means
    
    % Could use other methods too but bootsrapping appears to be most
    % adequate for non-normal distributions
    CI_temp=get_CI_bootstrap_tail(dat,25000); % get CI from bootstrap (25000 in paper)
    
    Low=CI_temp(1,:);
    Upp=CI_temp(2,:);
    
    % find significant FCCs
    isSignif=[Upp(m_vec_order)<0 | Low(m_vec_order)>0];
    
    % update variable names as we remove below tolerance variables
    dat_vars=variables(index);
    sort_vars=dat_vars(m_vec_order);
    sig_vars=sort_vars(isSignif);
    
    eval(['topVars.case',num2str(i),'=sig_vars(1:noTopVar);'])
end
% This is how we get our top candidatesusing the bootsrapping approach. 
varList=[topVars.case1;topVars.case2;topVars.case3;topVars.case4];
varList=unique(varList);

%% Bonferroni test for all cases
% Reorder variables according to magnitude in mean difference of case 1&2
% to help reading of the results
locVars=find_cell(varList,variables);

mC1=mean(case1(:,locVars));
mC2=mean(case2(:,locVars));
[~,idSortVarlist]=sort(abs(mC1-mC2),'descend');
varList=varList(idSortVarlist);

% lovate variables
locVars=find_cell(varList,variables);
% Perform ttest for each case with bonferroni correction
alpha=0.05;
pairCases=nchoosek(1:4,2); % combinations of FDP comparisons
ntest=size(pairCases,1)*numel(varList);
allCases=[];
figure;
allCIs=[];
myDiffs=[];
for pairNum=1:size(pairCases,1)
    tempCI=[];
    tempEst=[];
    for varNum=1:numel(varList)
        eval(['tempX=case',num2str(pairCases(pairNum,1)),'(:,locVars(varNum));'])
        eval(['tempY=case',num2str(pairCases(pairNum,2)),'(:,locVars(varNum));'])
        [H,P,CI,STATS]=ttest(tempX',tempY','alpha',alpha/ntest); % Note we divide by total number of test we compare overall (6 pairs * 15 vars)
        le_string=[varList{varNum},': C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))]; % Names of comparisons
        estimate=mean(tempX)-mean(tempY); %diff means
        % Collect data as we loop through the case comparisons
        allCases=[allCases;[{le_string},{estimate},{CI(1)},{CI(2)},{P}]];
        tempCI=[tempCI;CI];
        tempEst=[tempEst;estimate];
        myDiffs = [myDiffs tempX-tempY];
    end
    % Plot of comparisons pairwise
    subplot(3,2,pairNum)
    isSignif=0<tempCI(:,1)|0>tempCI(:,2);
    idSig=find(isSignif);
    idNotSig=find(~isSignif);
    errorbar(idSig,tempEst(idSig),tempCI(idSig,1)-tempEst(idSig)...
        ,tempEst(idSig)-tempCI(idSig,2),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',2);
    hold on
    errorbar(idNotSig,tempEst(idNotSig),tempCI(idNotSig,1)-tempEst(idNotSig)...
        ,tempEst(idNotSig)-tempCI(idNotSig,2),'db','MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    ylabel('Bonferroni CIs')
    title(['C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))])
    set(gca,'XTick',1:numel(varList))
    set(gca,'XTickLabel',varList)
    set(gca,'XTickLabelRotation',60)
    plot(xlim,[0 0],'--k')
    set(gca,'fontsize',22)
    set(gca,'fontweight','bold')
    grid on
    allCIs=[allCIs;tempCI];
end

%% Normal exact testing

% Get variance matrices and the means and std of the cases
Smat=zeros(4*numel(locVars)); %getting the variance matrix
vecMean=[];
vecSTD=[];
for caseNo=1:4
    eval(['tempData=case',num2str(caseNo),'(:,locVars);'])
    tempSmat=cov(tempData);
    idMat=(caseNo-1)*numel(locVars)+1:numel(locVars)*caseNo;
    Smat(idMat,idMat)=tempSmat;
    tempMean=mean(tempData);
    vecMean=[vecMean,tempMean];
    tempSTD=std(tempData);
    vecSTD=[vecSTD,tempSTD];
end
Smean=Smat/noSamples;

% Generate the covariance matrix for the case comparisons
pairCases=nchoosek(1:4,2);
ntest=size(pairCases,1)*numel(varList);
K=zeros(ntest,4*numel(varList));
countTests=0;
for pairNum=1:size(pairCases,1)
    for varNum=1:numel(varList)
        countTests=countTests+1;
        K(countTests,(pairCases(pairNum,1)-1)*numel(varList)+varNum)=1;
        K(countTests,(pairCases(pairNum,2)-1)*numel(varList)+varNum)=-1;
    end
end
Cont_est=K*vecMean';
Cont_var=K*Smean*transpose(K);
% Get correlation mat from the covariance that we just computed above
[Gamma,~] = corrcov(Cont_var);

nsimu=noSamples; % number of simulations
Zsimu = mvnrnd(zeros(ntest,1),Gamma,nsimu);% simulate Z nsimu times
absmax = max(abs(Zsimu)'); %compute max_j |Z_j|
q_alpha=quantile(absmax,1-alpha);

Low = Cont_est - q_alpha*sqrt(diag(Cont_var)); % We don't divide by n here as we did it earlier.
Upp = Cont_est + q_alpha*sqrt(diag(Cont_var));

figure
for pairNum=1:size(pairCases,1)
    tempEst=Cont_est((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum);
    tempCI=[Low((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum),...
        Upp((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum)];
    
    subplot(3,2,pairNum)
    isSignif=0<tempCI(:,1)|0>tempCI(:,2);
    idSig=find(isSignif);
    idNotSig=find(~isSignif);
    errorbar(idSig,tempEst(idSig),tempCI(idSig,1)-tempEst(idSig)...
        ,tempEst(idSig)-tempCI(idSig,2),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',2);
    hold on
    errorbar(idNotSig,tempEst(idNotSig),tempCI(idNotSig,1)-tempEst(idNotSig)...
        ,tempEst(idNotSig)-tempCI(idNotSig,2),'db','MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    ylabel('Exact Normal CIs')
    title(['C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))])
    set(gca,'XTick',1:numel(varList))
    set(gca,'XTickLabel',varList)
    set(gca,'XTickLabelRotation',60)
    plot(xlim,[0 0],'--k')
    set(gca,'fontsize',22)
    set(gca,'fontweight','bold')
    grid on
end

%% Bootstrap with tail balancing

% Generate the K matrix for matrix multiplications
pairCases=nchoosek(1:4,2);
ntest=size(pairCases,1)*numel(varList);
K=zeros(ntest,4*numel(varList));
countTests=0;
for pairNum=1:size(pairCases,1)
    for varNum=1:numel(varList)
        countTests=countTests+1;
        K(countTests,(pairCases(pairNum,1)-1)*numel(varList)+varNum)=1;
        K(countTests,(pairCases(pairNum,2)-1)*numel(varList)+varNum)=-1;
    end
end
vecMean=[];
vecSTD=[];
for caseNo=1:4
    eval(['tempData=case',num2str(caseNo),'(:,locVars);'])
    tempMean=mean(tempData);
    vecMean=[vecMean,tempMean];
    tempSTD=std(tempData);
    vecSTD=[vecSTD,tempSTD];
end

Boot_est=K*vecMean';
Boot_SD=sqrt(0.5*abs(K)*[vecSTD.^2]');

% Bootstrap samples
nR = 25000; 
t_stat = zeros(nR,ntest);
for i = 1:nR
    % Sampling with replacement
    index=randsample(noSamples,noSamples,true);
    tempDAT=[case1(index,locVars) case2(index,locVars) case3(index,locVars) case4(index,locVars)];
    tempDAT_m=mean(tempDAT);
    tempDAT_std=std(tempDAT);
    
    Samp_est=K*tempDAT_m';
    Samp_SD=sqrt(0.5*abs(K)*tempDAT_std.^2');
    t_stat(i,:)=(Samp_est-Boot_est)./(Samp_SD/sqrt(noSamples));
end

h_boot=[tiedrank(t_stat)-1]/[nR+1];
h_max=max(abs(h_boot-0.5)'); %take absolute here 

q_alpha=quantile(h_max,1-alpha);
beta = 1-2*q_alpha;

u_vec=quantile(t_stat,1-beta/2,1);
l_vec=quantile(t_stat,beta/2,1);

Low_tail = Boot_est - Boot_SD/sqrt(noSamples).*u_vec';
Upp_tail = Boot_est - Boot_SD/sqrt(noSamples).*l_vec';

figure
for pairNum=1:size(pairCases,1)
    tempEst=Boot_est((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum);
    tempCI=[Low_tail((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum),...
        Upp_tail((pairNum-1)*numel(locVars)+1:numel(locVars)*pairNum)];
    subplot(3,2,pairNum)
    errorbar(1:numel(tempEst),tempEst,tempEst-tempCI(:,1)...
        ,tempCI(:,2)-tempEst,'dk','MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5);
    isSignif=0<tempCI(:,1)|0>tempCI(:,2);
    idSig=find(isSignif);
    idNotSig=find(~isSignif);
    errorbar(idSig,tempEst(idSig),tempEst(idSig)-tempCI(idSig,1)...
        ,tempCI(idSig,2)-tempEst(idSig),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',2);
    hold on
    errorbar(idNotSig,tempEst(idNotSig),tempEst(idNotSig)-tempCI(idNotSig,1)...
        ,tempCI(idNotSig,2)-tempEst(idNotSig),'db','MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    ylabel('Bootstrap tail balance CIs')
    title(['C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))])
    set(gca,'XTick',1:numel(varList))
    set(gca,'XTickLabel',varList)
    set(gca,'XTickLabelRotation',60)
    hold on
    plot(xlim,[0 0],'--k')
    set(gca,'fontsize',22)
    set(gca,'fontweight','bold')
    grid on
end

%% Univariate test for all cases to show differences

% locate variables
locVars=find_cell(varList,variables);
% Perform ttest for each case with univariate method
alpha=0.05;
pairCases=nchoosek(1:4,2); % combinations of FDP comparisons
ntest=size(pairCases,1)*numel(varList);
allCases=[];
figure;
allCIs=[];
myDiffs=[];
for pairNum=1:size(pairCases,1)
    tempCI=[];
    tempEst=[];
    for varNum=1:numel(varList)
        eval(['tempX=case',num2str(pairCases(pairNum,1)),'(:,locVars(varNum));'])
        eval(['tempY=case',num2str(pairCases(pairNum,2)),'(:,locVars(varNum));'])
        le_string=[varList{varNum},': C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))]; % Names of comparisons
        estimate=mean(tempX-tempY); %mean of differences
        n=numel(tempX);
        CI = repmat(estimate,2,1) + [-1; 1]*tinv(1-alpha/2,n-1)*[std(tempX-tempY)/sqrt(n)];
        
        % Collect data as we loop through the case comparisons
        allCases=[allCases;[{le_string},{estimate},{CI(1)},{CI(2)},{P}]];
        tempCI=[tempCI;CI'];
        tempEst=[tempEst;estimate];
        myDiffs = [myDiffs tempX-tempY];
    end
    % Plot of comparisons pairwise
    subplot(3,2,pairNum)
    isSignif=0<tempCI(:,1)|0>tempCI(:,2);
    idSig=find(isSignif);
    idNotSig=find(~isSignif);
    errorbar(idSig,tempEst(idSig),tempCI(idSig,1)-tempEst(idSig)...
        ,tempEst(idSig)-tempCI(idSig,2),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',2);
    hold on
    errorbar(idNotSig,tempEst(idNotSig),tempCI(idNotSig,1)-tempEst(idNotSig)...
        ,tempEst(idNotSig)-tempCI(idNotSig,2),'db','MarkerFaceColor','w','MarkerSize',8,'LineWidth',2);
    ylabel('Univariate CIs')
    title(['C',num2str(pairCases(pairNum,1)),'-C',num2str(pairCases(pairNum,2))])
    set(gca,'XTick',1:numel(varList))
    set(gca,'XTickLabel',varList)
    set(gca,'XTickLabelRotation',60)
    plot(xlim,[0 0],'--k')
    set(gca,'fontsize',22)
    set(gca,'fontweight','bold')
    grid on
    allCIs=[allCIs;tempCI];
end

