% Script to show the ensemble estimation of Mnle.
% Used to make Figures 5 & S2.
clear;

% Define some constants.
Nc=1e4;
Ns=1e2;
m1=0.0;
m2=3.5;
b=1.0;
q=-1;
SMOOTHflag='none';

% Set the RNG seed in stone, so the plot is reproducible.
rng(2);

% Derive some variables.
T=0:Nc;
Mmax=linspace(0.6*m2,m2,Nc+1);
M=GR_MFD_Rand(m1,Mmax(2:end), log10(Nc),b, [1 Nc]);
Mlrg=OrderStatistic(M,Nc,'none'); Mlrg=[m1,Mlrg];
Mlrg_est=m1+log10(Nc)/b;

% Do the MLE-dM-fitting.
[m2_fit,llr]=M2fit(M,m1,b,1e2);

% Do the KS-dM-test (and regular KS-test).
Msamp=GR_MFD_Rand(m1,mean(m2), log10(Nc),b, [1 Nc]);
[~,KSp_fl]=kstest2(M,Msamp);
KSp_dm=KS_dM_test(M,m1,'resample');

% Propose possible three possible models of Mmax.
S(1).Mmax=Mmax; % True model.
S(2).Mmax=(Mmax(end)+0.0)*ones(size(Mmax));
S(3).Mmax=Inf*ones(size(Mmax));

% Get the ensemble weights.
W=EnsembleW(M,m1,S,b,SMOOTHflag);

% Get each Mmax model's estimate of the NLE's magnitude.
MnleE=zeros(size(Mlrg)); Wb=[];
for j=1:length(W)
    W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,b,q);
    MnleE=MnleE+(W(j).Mnle.*W(j).W);
    Wb=[Wb;W(j).W];
end

% Get the odds ratios.
ORm=Wb(1,end)/Wb(3,end);
ORg=Wb(2,end)/Wb(3,end);

% Report some values.
[KSp_fl KSp_dm]
[max(M) Mlrg_est mean(m2_fit) std(m2_fit)]
[ORm ORg]

% Define some colours I'd like to use.
C1='#345da7';
C2='#4c6ef5';
C3='#77aaff';
Rs=getMscale(M)/3;

% Plot.
figure(5); clf;
% MvT plot.
ax1=subplot(211);
scatter(T(2:end),M,Rs,'r','filled','HandleVisibility','off'); hold on;
plot(T,Mlrg ,'-r','HandleVisibility','off');
%plot(T,Mmax,'-b','HandleVisibility','off');
plot(T,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(T,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(T,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(T,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(T,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(T,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 max(T)]); ylim([m1 m2+0.5]);
legend('Location','southeast');
% WvT plot.
ax2=subplot(212);
bp=bar(T,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(S)
    plot(xlim(),[i i]/length(S),'-w');
end
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 max(T)]); ylim([0 1]);
linkaxes([ax1 ax2],'x');

% Plot again, but in log-scale this time.
figure(52); clf;
% MvT plot.
ax1=subplot(211);
scatter(T(2:end),M,Rs,'r','filled','HandleVisibility','off'); hold on;
plot(T,Mlrg ,'-r','HandleVisibility','off');
%plot(T,Mmax,'-b','HandleVisibility','off');
plot(T,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(T,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(T,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(T,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(T,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(T,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 max(T)]); ylim([m1 m2+0.5]);
set(gca, 'XScale', 'log');
legend('Location','southeast');
% WvT plot.
ax2=subplot(212);
bp=bar(T,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(S)
    plot(xlim(),[i i]/length(S),'-w');
end
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 max(T)]); ylim([0 1]);
set(gca, 'XScale', 'log');
linkaxes([ax1 ax2],'x');







%%%% SUBROUNTINE.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end