% Script to plot the case study data.
% See also script_CaseAnalysis.m for the analyses.
% Used to make Figure 8.
clear;

% Grab all of the cases to consider and load in.
%CaseList={'Basel','CB1','CB4','Paralana','Pohang','ESB10','PNR2','PeaceRiver','ParadoxValley','GGB','Prague','Pawnee','Youngstown'};
%D=PreProcData(CaseList);
load('CaseData_resamp.mat','D');

% Mmax model values.
Vx=10.^(1:0.1:9);
Mmax_Mx=Mmax_V(Vx,30,'McGarr');
Mmax_Gx=Mmax_V(Vx,1.5e9,'Galis');

% Make MvV plot.
figure(8); clf;
semilogx(Vx,Mmax_Mx,'-b'); hold on;
semilogx(Vx,Mmax_Gx,'--b');

% Preallocate.
ORm=zeros(size(D));
ORg=ORm;

% Loop over all of the case data.
for i=1:length(D)

    % Get the case of interest's catalogue.
    T=D(i).T;
    M=D(i).M;
    t=D(i).t;
    v=D(i).v;
    V=D(i).V;
    
    % Get some values.
    Nc=length(M);
    
    % Interpolate the injection volume.
    [t,I]=unique(t);
    V=V(I); v=v(I);
    Vi=interp1(t,V,T,'linear','extrap');
    
    % Mlrg sequence.
    Mlrg=OrderStatistic(M,Nc,'none');
    [Mlrg2,I]=unique(Mlrg);
    Nm=1:Nc;
    
    % Get the (log) odds ratio.
    ORm(i)=log10(D(i).W(1,end)/D(i).W(3,end));
    ORg(i)=log10(D(i).W(2,end)/D(i).W(3,end));
    
    % Add to the plot.
    plot(Vi,Mlrg,'-r');
    scatter(Vi(I),Mlrg2,[],mean([ORm(i),ORg(i)])*ones(size(Mlrg2)),'filled','MarkerEdgeColor','r');

    % Report some values.
    D(i).Case
    [geomean(D(i).KSp_fl) geomean(D(i).KSp_dm)]
    [max(D(i).M) mean(D(i).Mlrg_est) mean(D(i).m2) std(D(i).m2)]
    [10.^ORm(i) 10.^ORg(i)]
    
end
xlabel('Volume (m^3)'); ylabel('M_{LRG} Magnitude');
set(gca, 'XScale', 'log');
xlim([1e1 1e9]);
ylim([0 7]);
h = colorbar(); colormap(gca,R_colormap('Odds Ratio'));
clim([-1 +1]);
ylabel(h, 'Log_{10} Odds Ratio');



