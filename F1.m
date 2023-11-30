% Script to show order statistics applied to a synthetic MvT sequence.
% Used to make Figure 1.
clear;

% Define some constants.
Nc=1e4;
Ns=1e2;
m1=0.0;
m2=3.5;
b=1.0;

% Derive some more variables.
T=1:Nc;
Mmax=linspace(0.9*(m2-m1)+m1,m2,Nc);
M=GR_MFD_Rand(m1,Mmax, log10(Nc),b, [1 Nc]);
Rs=getMscale(M)/3;

% Get the sequence of largest events (and other order statistics).
Mlrg1=OrderStatistic(M,length(M)-0,'none'); %Mlrg1=[m1,Mlrg1];
Mlrg2=OrderStatistic(M,length(M)-1,'none'); %Mlrg2=[m1,Mlrg2];
Mlrg3=OrderStatistic(M,length(M)-2,'none'); %Mlrg3=[m1,Mlrg3];
Mlrg4=OrderStatistic(M,length(M)-3,'none'); %Mlrg4=[m1,Mlrg4];
dm_lrg=diff(unique(Mlrg1));

% Get a bunch more data.
for i=1:Ns
    Mi=GR_MFD_Rand(m1,Mmax, log10(Nc),b, [1 Nc]);
    Bi=diff(unique(OrderStatistic(Mi,length(M)-0,'none')));
    dm_lrg=[dm_lrg,Bi];
end

% Expected PDF
Mb=min(M):0.01:max(M);
[Mpdf,Mcdf,Msvf]=GR_MFD(Mb, m1,m2, log10(Nc),b,'normalized');

% Plot.
figure(1); clf;
subplot(2,2,[1 2]);
scatter(T,M,Rs,'r','filled'); hold on;
plot(T,Mmax,'-b');
plot(T,Mlrg1 ,'-r');
plot(T,Mlrg2,'-m');
plot(T,Mlrg3,'-g');
plot(T,Mlrg4,'-c');
xlabel('Time'); ylabel('Magnitude');
xlim([-max(T)/20 max(T)]);
subplot(223);
histogram(M,'Normalization','pdf'); hold on;
plot(Mb,Mpdf,'-k');
plot(sort(M),(1-(1:length(M))/length(M)),'-ob');
plot(Mb,Msvf,'-k');
xlabel('Magnitude'); ylabel('GR-MFD (PDF/CDF)');
set(gca,'YScale','log');
subplot(224);
histogram(dm_lrg,'Normalization','pdf'); hold on;
plot(Mb,Mpdf,'-k');
plot(sort(dm_lrg),(1-(1:length(dm_lrg))/length(dm_lrg)),'-ob');
plot(Mb,Msvf,'-k');
xlabel('Next Biggest Event Difference'); ylabel('GR-MFD (PDF/CDF)');
set(gca,'YScale','log');





%%%% SUBROUNTINE.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

