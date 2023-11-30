% Script to analyze the case data with Mmax tests (and to QC it).
% Used to make the data for Figure 8.
% Also used to make Figures S7-S32. 
clear;

% Define some constants.
Ns=5e1;
Rflag='resample';
Nr=1e2;
dMm=0.2;
dMb=0.1;
dMd=0.1;
q=-1;
SMOOTHflag='none';
%CaseList={'Basel','CB1','CB4','Paralana','Pohang','ESB10','PNR2','PeaceRiver','ParadoxValley','GGB','Prague','Pawnee','Youngstown'};
CaseList={'Paralana'};

% Get all of the case data.
D=PreProcData(CaseList);

% Loop over all of the k case data.
for k=1:length(D)
    
    % Get the case's catalogue.
    Lat=D(k).Lat'; Lon=D(k).Lon'; Dep=D(k).Dep';
    T=D(k).T'; M=D(k).M';
    T=datetime(T,'ConvertFrom','datenum');
    m1b=D(k).Mc;
    m1k=D(k).Mk;
    
    % Get the case's injection information.
    t=D(k).t'; t=datetime(t,'ConvertFrom','datenum');
    v=D(k).v';
    V=D(k).V';
    
    % Loop over the bootstrap trials.
    for i=1:Ns
        i
        
        % Apply a perturbation.
        Mi=M+(dMd*rand(size(M))-dMd/2);
        m1bi=m1b+((dMd/2)*rand(size(m1b))-dMd/4);
        m1ki=m1k+((dMd/2)*rand(size(m1k))-dMd/4);
        
        % Truncate on the magnitude of completeness.
        Imb=(Mi>=m1bi);
        Imk=(Mi>=m1ki);
        %ID=ID(Im); Lat=Lat(Im); Lon=Lon(Im); Dep=Dep(Im); T=T(Im); M=M(Im);
        
        % Derive some variables.
        Ncb=length(Mi(Imb));
        Nck=length(Mi(Imk));
        
        % Fit the GR-MFD b-value.
        [bi,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(Mi,m1bi,dMb);
        b(i)=bi;
        
        % Get the expected Mlrg value.
        Mlrg_est(i)=m1bi+log10(Ncb)/bi;
        
        % Do the MLE-dM-fitting.
        [m2_fit,llr]=M2fit(Mi(Imb),m1bi,bi,Nr);
        m2_avg(i)=mean(m2_fit);
        m2_std(i)=std(m2_fit);
        m2_p90(i)=prctile(m2_fit,90);
        
        % Do the KS-dM-test (and regular KS-test).
        Msamp=GR_MFD_Rand(m1ki,m2_avg(i), log10(Nck),bi, [1 Nck]);
        [~,p_full]=kstest2(Mi(Imk),Msamp);
        KSp_fl(i)=p_full;
        KSp_dm(i)=KS_dM_test(Mi(Imk),m1ki,Rflag);
        
    end
    
    % Save data into the output structure.
    D(k).b=b;
    D(k).Mlrg_est=Mlrg_est;
    D(k).m2=m2_fit;
    D(k).KSp_fl=KSp_fl;
    D(k).KSp_dm=KSp_dm;
    
    % Truncate on the magnitude of completeness.
    Imb=(M>=m1b);
    Ncb=length(M(Imb));
    
    % Get the sequence of largest events (and indicies to when they occur).
    Mlrg=OrderStatistic(M(Imb),length(M(Imb))-0,'none'); Mlrg=[m1b,Mlrg];
    [~,Id]=unique(Mlrg);
    Tlrg=[T(find(Imb,1,'first'))-1/24/3600,T(Imb)];
    
    % Interpolate
    [t,I]=unique(t);
    V=V(I); v=v(I);
    Vi=interp1(t,V,Tlrg,'linear','extrap');

    % Ad hoc fix for EQs before injection start with CB1 and Pohang.
    Vi(Vi==0)=min(Vi(Vi>0));
    
    % Fit the McGarr Mmax model.
    Afxn = @(G) min(Mmax_V(Vi,G,'McGarr')-Mlrg)-dMm;
    G=fzero(Afxn,[1e-2 1e5]);
    Mmax_M=Mmax_V(Vi,G,'McGarr');

    % Fit the Galis Mmax model.
    Afxn = @(g) min(Mmax_V(Vi,g,'Galis')-Mlrg)-dMm;
    g=fzero(Afxn,[1e-2 1e15]);
    Mmax_G=Mmax_V(Vi,g,'Galis');
    
    % Regular Mmax models.
    %Mmax_M=Mmax_V(Vi,30,'McGarr');
    %Mmax_G=Mmax_V(Vi,1.5e9,'Galis');
    
    % Propose possible three possible models of Mmax.
    S(1).Mmax=Mmax_M;
    S(2).Mmax=Mmax_G;
    S(3).Mmax=(Inf)*ones(size(Mlrg));
    
    % Get the ensemble weights.
    W=EnsembleW(M(Imb),m1b,S,mean(bi),SMOOTHflag);
    
    % Get each Mmax model's estimate of the NLE's magnitude.
    MnleE=zeros(size(Mlrg)); Wb=[];
    for j=1:length(W)
        W(j).Mnle=NLE_M(Mlrg,W(j).Mmax,mean(bi),q);
        MnleE=MnleE+(W(j).Mnle.*W(j).W);
        Wb=[Wb;W(j).W];
    end
    
    % Save some more data into the output structure.
    D(k).G=G;
    D(k).g=g;
    D(k).W=Wb;
    
    % Get the GR-MFD plotting values.
    po=[-mean(b),a];
    Mgr_fit=[m1b, max(M)];
    Ngr_fit=10.^polyval(po,Mgr_fit);
    
    % Report the case name
    D(k).Case
    
    % Report the KS-dM-test p-value.
    geomean(KSp_dm)
    
    % Report the MLE-fitted Mmax (and related) values.
    mean(Mlrg_est)
    mean(m2_avg)
    mean(m2_std)
    
    % Report the model weights and the relative odds of each model.
    Wb(:,end)
    [Wb(:,end)/Wb(1,end),Wb(:,end)/Wb(2,end),Wb(:,end)/Wb(3,end)]
    
end

% Save data.
save('CaseData_temp.mat','D');

% Define some colours I'd like to use.
C1='#345da7';
C2='#4c6ef5';
C3='#77aaff';
GREY=[0.85,0.85,0.85];
Rs=getMscale(M)/3;

% Plot.
figure(1); clf;
% MvT plot.
ax1=subplot(211);
scatter(1:Ncb,M(Imb),Rs(Imb),'r','filled','HandleVisibility','off'); hold on;
%scatter(T(~Im),M(~Im),Rs(~Im),'m','filled','HandleVisibility','off');
plot(0:Ncb,Mlrg ,'-r','HandleVisibility','off');
plot(0:Ncb,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(0:Ncb,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(0:Ncb,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(0:Ncb,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(0:Ncb,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(0:Ncb,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([0 Ncb]); %ylim([m1 4.2]);
legend('Location','southeast');
% WvT plot.
ax2=subplot(212);
bp=bar(0:Ncb,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
for i=1:length(S)
    plot(xlim(),[i i]/length(S),'-w');
end
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([0 Ncb]); ylim([0 1]);
linkaxes([ax1 ax2],'x');

% Plot.
figure(58); clf;
% MvT plot.
ax1=subplot(211);
scatter(1:Ncb,M(Imb),Rs(Imb),'r','filled','HandleVisibility','off'); hold on;
plot(0:Ncb,Mlrg ,'-r','HandleVisibility','off');
plot(0:Ncb,W(1).Mmax,'-b', 'DisplayName','M_{MAX} Model 1');
plot(0:Ncb,W(2).Mmax,'--b','DisplayName','M_{MAX} Model 2');
plot(0:Ncb,W(3).Mmax,':b', 'DisplayName','M_{MAX} Model 3');
plot(0:Ncb,W(1).Mnle,'-g', 'DisplayName','M_{NLE} Model 1');
plot(0:Ncb,W(2).Mnle,'--g','DisplayName','M_{NLE} Model 2');
plot(0:Ncb,W(3).Mnle,':g', 'DisplayName','M_{NLE} Model 3');
xlabel('Chronological Event Number'); ylabel('Magnitude');
xlim([1 Ncb]);
set(gca, 'XScale', 'log');
legend('Location','northwest');
% WvT plot.
ax2=subplot(212);
bp=bar(0:Ncb,Wb,'stacked','LineStyle','none','BarWidth',1); hold on;
bp(1).FaceColor=C1; bp(2).FaceColor=C2; bp(3).FaceColor=C3;
for i=1:length(S)
    plot([1 Ncb],[i i]/length(S),'-w');
end
xlabel('Chronological Event Number'); ylabel('Model Weights');
xlim([1 Ncb]); ylim([0 1]);
set(gca, 'XScale', 'log');
linkaxes([ax1 ax2],'x');

% Plot catalogue filtering info.
figure(57); clf;
% Map.
ax1=subplot(3,3,[1 2 4 5]);
if(~isempty(Lat))
    scatter(Lon(~Imb),Lat(~Imb),Rs(~Imb),'m','filled'); hold on;
    scatter(Lon(Imb), Lat(Imb), Rs(Imb), 'r','filled');
end
xlabel('Longitude'); ylabel('Latitude');
% Depth.
ax2=subplot(3,3,[3 6]);
if(~isempty(Lat))
    scatter(Dep(~Imb),Lat(~Imb),Rs(~Imb),'m','filled'); hold on;
    scatter(Dep(Imb), Lat(Imb), Rs(Imb), 'r','filled');
end
xlabel('Depth (km)'); ylabel('Latitude');
xlim([0 10]);
% GR-FMD.
ax3=subplot(3,3,[7 8 9]);
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlim([min(Mgr)-dMd/2 max(Mgr)+dMd/2]); ylim([0.7 1.3*max(Ngr)]);
plot(m1b*[1 1],ylim,'--k');
xlabel('Magnitude (M_L)'); ylabel('Count');






%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(Mw)
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end