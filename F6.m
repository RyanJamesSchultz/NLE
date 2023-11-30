% Script to apply Mmax tests to the Groningen data.
% Used to make Figures 6, S3, & S4.
clear;

% Define some values.
Ns=5e1;
m1=1.5;
dM=0.1;
Rflag='resample';
Nr=1e2;
q=-1;
SMOOTHflag='none';

% Set the RNG seed in stone, so the plot is reproducible.
rng(6);

% Get the Groningen catalogue.
load('data/GRONINGEN_cat.mat','Catalog'); 
ID=[Catalog(1).val]'; Lat=[Catalog(3).val]'; Lon=[Catalog(4).val]'; Dep=[Catalog(5).val]';
T=[Catalog(2).val]'; T=datetime(T,'ConvertFrom','datenum');
M=[Catalog(6).val]';

% Fitler the catalogue spatiotemporally.
[ID,Lon,Lat,Dep,T,M,~]=filtCat(ID,Lon,Lat,Dep,T,M);

% Loop over the bootstrap trials.
for i=1:Ns
    i
    
    % Apply a perturbation.
    Mi=M+(dM*rand(size(M))-dM/2);
    m1i=m1+((dM/2)*rand(size(m1))-dM/4);
    
    % Truncate on the magnitude of completeness.
    Im=(Mi>=m1i);
    %ID=ID(Im); Lat=Lat(Im); Lon=Lon(Im); Dep=Dep(Im); T=T(Im); M=M(Im);
    
    % Derive some variables.
    Nc=length(Mi(Im));
    
    % Fit the GR-MFD b-value.
    [bi,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(Mi,m1i,dM);
    b(i)=bi;
    
    % Get the expected Mlrg value.
    Mlrg_est(i)=m1i+log10(Nc)/bi;
    
    % Do the MLE-dM-fitting.
    [m2_fit,llr]=M2fit(Mi(Im),m1i,bi,Nr);
    m2_avg(i)=mean(m2_fit);
    m2_std(i)=std(m2_fit);
    m2_p90(i)=prctile(m2_fit,90);
    
    % Do the KS-dM-test.
    Msamp=GR_MFD_Rand(m1i,m2_avg(i), log10(Nc),bi, [1 Nc]);
    [~,p_full]=kstest2(Mi(Im),Msamp);
    KSp_fl(i)=p_full;
    KSp_dm(i)=KS_dM_test(Mi(Im),m1i,Rflag);
    
end

% Truncate on the magnitude of completeness.
Im=(M>=m1);
Nc=length(M(Im));

% Get the sequence of largest events (and indicies to when they occur).
Mlrg=OrderStatistic(M(Im),length(M(Im))-0,'none'); Mlrg=[m1,Mlrg];
[~,Id]=unique(Mlrg);
Tlrg=[T(find(Im,1,'first'))-1/24/3600,T(Im)];

% Propose possible three possible models of Mmax.
S(1).Mmax=(mean(m2_avg))*ones(size(Mlrg));
S(2).Mmax=(mean(m2_avg)+mean(m2_std))*ones(size(Mlrg));
S(3).Mmax=(Inf)*ones(size(Mlrg));

% Get the ensemble weights.
W=EnsembleW(M(Im),m1,S,mean(bi),SMOOTHflag);

% Get each Mmax model's estimate of the NLE's magnitude.
MnleE=zeros(size(Mlrg)); Wb=[];
for j=1:length(W)
    W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,mean(bi),q);
    MnleE=MnleE+(W(j).Mnle.*W(j).W);
    Wb=[Wb;W(j).W];
end

% Get the GR-MFD plotting values.
po=[-mean(b),a];
Mgr_fit=[m1, max(M)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Report the KS-dM-test p-value.
geomean(KSp_dm)

% Report the MLE-fitted Mmax (and related) values.
mean(Mlrg_est)
mean(m2_avg)
mean(m2_std)

% Report the model weights and the relative odds of each model. 
Wb(:,end)
[Wb(:,end)/Wb(1,end),Wb(:,end)/Wb(2,end),Wb(:,end)/Wb(3,end)]




% Define some colours I'd like to use.
C1='#345da7';
C2='#4c6ef5';
C3='#77aaff';
GREY=[0.85,0.85,0.85];
Rs=getMscale(M)/3;

% Plot.
figure(6); clf;
% MvT plot.
ax1=subplot(211);
scatter(1:Nc,M(Im),Rs(Im),'r','filled','HandleVisibility','off'); hold on;
%scatter(T(~Im),M(~Im),Rs(~Im),'m','filled','HandleVisibility','off');
plot(0:Nc,Mlrg ,'-r','HandleVisibility','off');
plot(0:Nc,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(0:Nc,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(0:Nc,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(0:Nc,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(0:Nc,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(0:Nc,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([0 Nc]); ylim([m1 4.2]);
legend('Location','southeast');
% WvT plot.
ax2=subplot(212);
bp=bar(0:Nc,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(S)
    plot(xlim(),[i i]/length(S),'-w');
end
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([0 Nc]); ylim([0 1]);
linkaxes([ax1 ax2],'x');

% Plot.
figure(54); clf;
% MvT plot.
ax1=subplot(211);
scatter(1:Nc,M(Im),Rs(Im),'r','filled','HandleVisibility','off'); hold on;
plot(0:Nc,Mlrg ,'-r','HandleVisibility','off');
plot(0:Nc,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(0:Nc,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(0:Nc,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(0:Nc,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(0:Nc,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(0:Nc,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 Nc]);
set(gca, 'XScale', 'log');
legend('Location','southeast');
% WvT plot.
ax2=subplot(212);
bp=bar(0:Nc,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
for i=1:length(S)
    plot([1 Nc],[i i]/length(S),'-w');
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 Nc]); ylim([0 1]);
set(gca, 'XScale', 'log');
linkaxes([ax1 ax2],'x');

% Plot catalogue filtering info.
figure(53); clf;
% Map.
ax1=subplot(3,3,[1 2 4 5]);
scatter(Lon(~Im),Lat(~Im),Rs(~Im),'m','filled'); hold on;
scatter(Lon(Im), Lat(Im), Rs(Im), 'r','filled');
xlabel('Longitude'); ylabel('Latitude');
% Depth.
ax2=subplot(3,3,[3 6]);
scatter(Dep(~Im),Lat(~Im),Rs(~Im),'m','filled'); hold on;
scatter(Dep(Im), Lat(Im), Rs(Im), 'r','filled');
xlabel('Depth (km)'); ylabel('Latitude');
xlim([0 4]);
% GR-FMD.
ax3=subplot(3,3,[7 8 9]);
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);
plot(m1*[1 1],ylim,'--k');
xlabel('Magnitude (M_L)'); ylabel('Count');






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

% Filter the catalogue spatiotemporally.
function [id,x,y,z,t,m,I,m1]=filtCat(id,x,y,z,t,m)
  %I=0;
  
  % Define filter details.
  z_L=[0 10];
  x_L=[ 6.5  7.0  7.0  6.5];
  y_L=[53.5 53.5 53.1 53.1];
  T_L=[datetime(1980,01,01) datetime(2023,01,01)];

  % Filter laterally, temporally and in depth.
  Is=inpolygon(x,y,x_L,y_L);
  It=(t>=min(T_L))&(t<=max(T_L));
  Id=(z>=min(z_L))&(z<=max(z_L));
  
  % Apply the filter.
  I=Is&Id&It;
  id=id(I);
  x=x(I);
  y=y(I);
  z=z(I);
  t=t(I);
  m=m(I);
end