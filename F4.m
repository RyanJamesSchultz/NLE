% Script to show that dMlrg is more sensitive to Mmax than the full catalogue.
% Used to make Figure 4.
clear;

% Define some variables.
Nc=1e3;
Ns=5e2;
Rflag='resample';
Mbig=3.0;
dMbig=-1.0;
b=1.0;

% Compute some dependent values.
Mmax=Mbig+dMbig;
m1=Mbig-log10(Nc)/b;
m2=Mmax;
a=log10(Nc);
T=1:Nc;

% Preallocate.
Pfull=zeros([1 Ns]);
Pmlrg=Pfull; dP=Pfull;

% Loop.
for i=1:Ns
    i
    
    % Get the catalogue.
    M=GR_MFD_Rand(m1,m2, a,b, [1 Nc]);
    
    % Do the KS-tests.
    Msamp=GR_MFD_Rand(m1,Inf, a,b, [1 Nc]);
    [~,p_full]=kstest2(M   ,Msamp);
    [p_mlrg]=KS_dM_test(M,m1,Rflag);
    p_mlrg=mean(p_mlrg);
    
    % Save results in vectors.
    Pfull(i)=log10(p_full);
    Pmlrg(i)=log10(p_mlrg);
    dP(i)=log10(p_full)-log10(p_mlrg);
end

% Fraction better than the control test.
sum(dP>0)/length(dP)

% Plot.
GREY=[0.85,0.85,0.85];
figure(4); clf;
subplot(121);
histogram(Pfull,'FaceColor',GREY); hold on;
histogram(Pmlrg);
plot(log10(0.05)*[1 1],ylim(),':k');
xlabel('log_{10} p-value'); ylabel('Counts');
subplot(122);
histogram(dP); hold on;
plot(mean(dP)*[1 1],ylim(),'-b');
plot(0*[1 1],ylim(),':k');
xlabel('log_{10} p-value difference'); ylabel('Counts');


