% Numerical 'proof' that the PDF/CDF for a GR-MFD on magnitudes is the same as the PDF/CDF for dM_lrg.
% Used to make Figure 2.
clear;

% Define some constants.
N=1e4;
Ns=1e3;
m1=0.0;
m2=2.0;
b=1.0;

% Derive some more variables.
a=log10(N);
T=1:N;

% Preallocate.
MlrgU=[];   MlrgB=[];   
dM_lrgU=[]; dM_lrgB=[]; 

% Loop.
for i=1:Ns

    % Draw a catalogue.
    Mb=GR_MFD_Rand(m1,m2 , a,b, [1 N]);
    Mu=GR_MFD_Rand(m1,Inf, a,b, [1 N]);
    
    % Get the biggest event sequence (for the bounded catalogue).
    Mlrg1_i=OrderStatistic(Mb,length(Mb)-0,'unique'); Mlrg1_i=[m1,Mlrg1_i];
    MlrgB=[MlrgB, Mlrg1_i];
    dM_lrgB=[dM_lrgB, diff(Mlrg1_i)];

    % Get the biggest event sequence (for the unbounded catalogue).
    Mlrg1_i=OrderStatistic(Mu,length(Mu)-0,'unique'); Mlrg1_i=[m1,Mlrg1_i];
    MlrgU=[MlrgU, Mlrg1_i];
    dM_lrgU=[dM_lrgU, diff(Mlrg1_i)];
    
end

% Test the sample distributions for similarity.
[~,Pu ]=kstest2(dM_lrgU ,Mu-m1) % Indistinguishable.
[~,Pb ]=kstest2(dM_lrgB ,Mb-m1) % Distinguishable.

% Expected PDF
mrU=min(Mb):0.1:max([Mu dM_lrgU]);
mrB=min(Mb):0.1:m2;
[MpdfU,McdfU,MsvfU]=GR_MFD(mrU, m1,Inf, a,b,'normalized');
[MpdfB,McdfB,MsvfB]=GR_MFD(mrB, m1, m2, a,b,'normalized');

% Plot.
figure(2); clf;
% Regular unbounded GR-MFD.
subplot(2,2,1);
histogram(Mu-m1,mrU-m1,'Normalization','pdf'); hold on;
plot(mrU-m1,MpdfU,'-k');
plot(mrU-m1,100*MsvfU,'-k','LineWidth',2);
plot(sort(Mu-m1),100*(1-(1:length(Mu))/length(Mu)),'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]);
YL=ylim();
% Unbounded dMlrg.
subplot(2,2,2);
histogram(dM_lrgU,mrU-m1,'Normalization','pdf'); hold on;
plot(mrU-m1,MpdfU,'-k');
plot(mrU-m1,100*MsvfU,'-k','LineWidth',2);
plot(sort(dM_lrgU),100*(1-(1:length(dM_lrgU))/length(dM_lrgU)),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);
% Regular bounded GR-MFD.
subplot(2,2,3);
histogram(Mb-m1,mrU-m1,'Normalization','pdf'); hold on;
plot(mrB-m1,MpdfB,'-k');
plot(mrB-m1,100*MsvfB,'-k','LineWidth',2);
plot(sort(Mb-m1),100*(1-(1:length(Mb))/length(Mb)),'-b','LineWidth',2);
xlabel('M-m_1'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
plot((m2-m1)*[1 1],ylim(),':k');
xlim([0 4]); ylim(YL);
% Bounded dMlrg.
subplot(2,2,4);
histogram(dM_lrgB,mrU-m1,'Normalization','pdf'); hold on;
plot(mrB-m1,MpdfB,'-k');
plot(mrB-m1,100*MsvfB,'-k','LineWidth',2);
plot(sort(dM_lrgB),100*(1-(1:length(dM_lrgB))/length(dM_lrgB)),'-b','LineWidth',2);
xlabel('\DeltaM_{LRG}'); ylabel('GR-MFD (PDF/1-CDF)');
set(gca,'YScale','log');
xlim([0 4]); ylim(YL);





