% Script that shows the likelihood ratio function for m2 fitting.
% Used to make Figure S1.
clear;

% Define some variables.
Nc=1e3;
Mbig=4.5;
dMbig=-1.0;
b=1.0;
Ns=1e2;

% Compute some dependent values.
Mmax=Mbig+dMbig;
m1=Mbig-log10(Nc)/b;
m2=Mmax;
a=log10(Nc);
T=1:Nc;

% Get a catalogue.
M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);
[Mlrg]=OrderStatistic(M,Nc,'unique'); Mlrg=[m1,Mlrg]; dMlrg=diff(Mlrg); Mlrg(end)=[];
M2=max(M):0.01:max(M)+2;

% Get the log-likelihood ratio function.
for j=1:length(M2)
    LLR(j)=+2*GR_MFD_LL(dMlrg,0,M2(j)-Mlrg,a,b)-2*GR_MFD_LL(M,m1,M2(j),a,b)-2*GR_MFD_LL(dMlrg,0,Inf,a,b)+2*GR_MFD_LL(M,m1,Inf,a,b);
end

% Preallocate.
Mlrg_t=Mlrg;
dMlrg_t=dMlrg;

% Loop.
for i=1:Ns
    i
    % Get a newly reshuffled catalogue.
    In=randperm(Nc);
    Mi=M(In);
    
    % Get the Mlrg sequence.
    [Mlrg_i]=OrderStatistic(Mi,Nc,'unique'); Mlrg_i=[m1,Mlrg_i]; dMlrg_i=diff(Mlrg_i); Mlrg_i(end)=[];
    
    % Check if this is the same as fitting the data. 
    [v1,v2]=M2fit(Mi,m1,b,0);
    M2_fit(i)=v1;
    LLR_fit(i)=v2;
    
    % Loop over all of the possible Mmax values.
    for j=1:length(M2)
        
        % Get log-likelihood ratio scores.
        llr(j,i)=+2*GR_MFD_LL(dMlrg_i,0,M2(j)-Mlrg_i,a,b)-2*GR_MFD_LL(Mi,m1,M2(j),a,b)-2*GR_MFD_LL(dMlrg_i,0,Inf,a,b)+2*GR_MFD_LL(M,m1,Inf,a,b);
    end
    
end

% Get the average estimates.
[R,Im]=min(llr,[],1);
Mfit=M2(Im);
Mfit_avg=mean(Mfit);

% Get the min LLR estimates.
for i=1:length(Im)
    LLfit(i)=llr(Im(i),i);
end

% Report values.
Mmax
mean(Mfit)
median(Mfit)

% Plot.
figure(51); clf;
% Plot the resampled log-likelihood ratio functions.
subplot(211);
plot(M2,llr,'-r'); hold on;
plot(M2,mean(llr,2),'-k','LineWidth',2);
plot(Mfit_avg*[1 1],ylim(),'--r','LineWidth',2);
plot(Mmax*[1 1],ylim(),'-k');
xlabel('Possible M_{MAX}'); ylabel('Log-Likelihood Ratio');
xlim([min(M2) max(M2)]);
% Plot the histogram of fitted Mmax values.
subplot(212);
histogram(Mfit,max(M):0.02:max(M)+2,'FaceColor','r'); hold on;
plot(Mfit_avg*[1 1],ylim(),'--r','LineWidth',2);
plot(Mmax*[1 1],ylim(),'-k');
xlabel('Possible M_{MAX}'); ylabel('Counts');
xlim([min(M2) max(M2)]);

% Checking the FS1 plots are consistent with the M2fit routine.
figure(2); clf;
% scatter plots.
subplot(221);
plot(Mfit,M2_fit,'o'); hold on;
plot([min(Mfit) max(Mfit)],[min(Mfit) max(Mfit)],'-k');
xlabel('Plot fit M2'); ylabel('Algo fit M2');
subplot(222);
plot(LLfit,LLR_fit,'o'); hold on;
plot([min(LLfit) max(LLfit)],[min(LLfit) max(LLfit)],'-k');
xlabel('Plot fit LLR'); ylabel('Algo fit LLR');
% residual histograms.
subplot(223);
histogram(Mfit-M2_fit);
xlabel('Mag difference'); ylabel('Counts');
subplot(224);
histogram(LLfit-LLR_fit);
xlabel('LLR difference'); ylabel('Counts');


