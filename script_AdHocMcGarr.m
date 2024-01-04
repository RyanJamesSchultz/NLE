% Script to test how well the ensemble weighting will discern McGarr/Galis 
% models if the relationships are fitted after the data has been drawn.  
% Also used to do some input sensitivity tests to examine for output result 
% biasing.
% 
% Used to derive results for Supplementary Text S4.
clear;

% Define some constants.
SI=-2.0;
b=1.0;
db=normrnd(0.0,0.05,1);
Rf=0.1;
dMd=0.1;
m1=0.0;
Vt=1e6;
dM=0.2;
q=-1;
Nr=1e2;
Rflag='resample';
SMOOTHflag='none';

% Set the RNG seed in stone, for reproducible tests.
%rng(2);

% Derive some catalogue-related variables.
Nc=Vt*10^(SI-b*m1);
M=GR_MFD_Rand(m1,Inf,log10(Nc),b,[1 Nc]);
Vi=linspace(0,Vt,Nc+1); Vi(1)=[];

% Perturb the b-value or round the catalogue magnitudes.
b=b+db;
if(Rf~=0)
    M=round(M/Rf)*Rf;
end

% Undo the rounding of magnitudes by dithering.
M=M+(dMd*rand(size(M))-dMd/2);
%M=M+(GR_MFD_Rand(0,dMd,0,b,size(M))-dMd/2);
m1=min(M);

% Get the sequence of largest events (and indicies to when they occur).
Mlrg=OrderStatistic(M,length(M)-0,'none'); Mlrg=[m1,Mlrg];
[~,Id]=unique(Mlrg);

% Do the KS-dM-test.
KSp=KS_dM_test(M,m1,Rflag);

% Do the MLE-dM-fitting.
[m2_fit,llr]=M2fit(M,m1,b,Nr);

% Fit the McGarr Mmax model.
Afxn = @(G) min(Mmax_V(Vi,G,'McGarr')-M)-dM;
G=fzero(Afxn,[1e-20 1e5]);
%G=30;
Mmax_M=Mmax_V(Vi,G,'McGarr');

% Fit the Galis Mmax model.
Afxn = @(g) min(Mmax_V(Vi,g,'Galis')-M)-dM;
g=fzero(Afxn,[1e-20 1e12]);
%g=1.5e9;
Mmax_G=Mmax_V(Vi,g,'Galis');

% Propose three Mmax models.
S(1).Mmax=[Mmax_M(1),Mmax_M];
S(2).Mmax=[Mmax_M(1),Mmax_G];
S(3).Mmax=(Inf)*ones(size(Mlrg));

% Get the ensemble weights.
W=EnsembleW(M,m1,S,b,SMOOTHflag);

% Get each Mmax model's estimate of the NLE's magnitude.
MnleE=zeros(size(Mlrg)); Wb=[];
for j=1:length(W)
    W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,b,q);
    MnleE=MnleE+(W(j).Mnle.*W(j).W);
    Wb=[Wb;W(j).W];
end

% Define some colours I'd like to use.
C1='#345da7';
C2='#4c6ef5';
C3='#77aaff';
GREY=[0.85,0.85,0.85];
Rs=getMscale(M)/3;

% Report values.
G
g
KSp
mean(m2_fit)
std(m2_fit)
Wb(:,end)
[Wb(:,end)/Wb(1,end),Wb(:,end)/Wb(2,end),Wb(:,end)/Wb(3,end)]

% Plot.
figure(1); clf;
% MvT plot.
ax1=subplot(211);
scatter(Vi,M,Rs,'r','filled','HandleVisibility','off'); hold on;
plot([0,Vi],Mlrg ,'-r','HandleVisibility','off');
plot([0,Vi],W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot([0,Vi],W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot([0,Vi],W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot([0,Vi],W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot([0,Vi],W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot([0,Vi],W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Injected Volume (m^3)'); ylabel('Magnitude');
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvT plot.
ax2=subplot(212);
bp=bar([0,Vi],Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
set(gca, 'XScale', 'log');
for i=1:length(S)
    plot([min(Vi) max(Vi)],[i i]/length(S),'-w');
end
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([min(Vi) max(Vi)]); ylim([0 1]);
linkaxes([ax1 ax2],'x');

% Fit the GR-MFD b-value.
%[bf,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(M,m1,Rf);




%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end
