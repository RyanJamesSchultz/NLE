% Loads the test bench data to make Tables 3 & 4.
% See also script_TestBench.m
clear;

% Load in the data
load('TestBench_resamp.mat','S');
Nc   =reshape([S.Nc]   ,size(S));
dMbig=reshape([S.dMbig],size(S));

% Get some information.
Ns=length(S(1,1).Mmax);
Nh=round(1*sqrt(Ns));

% Grab values of interest.
for i=1:size(S,1)
    for j=1:size(S,2)
        Pfull(i,j)=mean(S(i,j).Pfull);
        Pmlrg(i,j)=mean(S(i,j).Pmlrg);
        dPm(i,j)=mean(S(i,j).Pfull-S(i,j).Pmlrg);
        dM(i,j)=median(S(i,j).Mmax)-(4+dMbig(i,j));
        dMv(i,j)=std(S(i,j).Mmax);
        dMe(i,j)=mean(S(i,j).Mmax_err);
    end
end

% Select table entries to plot in greater detail.
i=3;
j=8;
Nc(i,j)
dMbig(i,j)

% Plot KS-test stuff.
dP=S(i,j).Pfull-S(i,j).Pmlrg;
GREY=[0.85,0.85,0.85];
figure(1); clf;
subplot(121);
histogram(S(i,j).Pfull,Nh,'FaceColor',GREY); hold on;
histogram(S(i,j).Pmlrg,Nh);
plot(log10(0.05)*[1 1],ylim(),':k');
xlabel('log_{10} p-value'); ylabel('Counts');
subplot(122);
histogram(dP); hold on;
plot(mean(dP)*[1 1],ylim(),'-b');
plot(0*[1 1],ylim(),':k');
xlabel('log_{10} p-value difference'); ylabel('Counts');

% Plot MLE stuff.
figure(2); clf;
histogram(S(i,j).Mmax,Nh,'FaceColor','r'); hold on;
plot((4+dMbig(i,j))*[1 1],ylim(),':k');
plot(mean([S(i,j).Mmax])*[1 1],ylim(),'-r');
xlabel('M_{MAX}'); ylabel('Counts');

% Table 3 details.
Pmlrg
Pfull

% Table 4 details.
dM
