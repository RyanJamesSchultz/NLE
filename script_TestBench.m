% Test bench to explore how well KS(dM)-test and m2-fitting works.
% Used to make the data for Tables 3 & 4.
clear;

% Define some variables.
Ns=5e2;
Nc=[1e2 1e3 1e4 1e5 1e6];
dMbig=[-2.0 -1.5 -1.0 -0.5 0.0 +0.5 +1.0 Inf];
Mbig=4.0;
b=1.0;
Rflag='resample';
Nr=1e2;

% Preallocate.
S=struct('Nc',0,'dMbig',0,'Pfull',zeros([1 Ns]),'Pmlrg',zeros([1 Ns]),'Mmax',zeros([1 Ns]),'Mmax_err',zeros([1 Ns]));

% Loop over parameter space.
for j=1:length(Nc)
    for k=1:length(dMbig)
        
        % Save some important details.
        S(j,k).Nc=Nc(j);
        S(j,k).dMbig=dMbig(k);
        
        % Compute some dependent values.
        Mmax=Mbig+dMbig(k);
        m1=Mbig-log10(Nc(j))/b;
        m2=Mmax;
        a=log10(Nc(j));
        
        % Loop over all of the catalogue trials.
        for i=1:Ns
            [i,k,j]
            
            % Get the catalogue.
            M=GR_MFD_Rand(m1,m2, a,b, [1 Nc(j)]);
            Msamp=GR_MFD_Rand(m1,Inf, a,b, [1 1e6]);
            
            % Do the KS-tests.
            [~,p_full]=kstest2(M,Msamp);
            p_mlrg=KS_dM_test(M,m1,Rflag);
            
            % Fit m2.
            [m2_fit,llr]=M2fit(M,m1,b,Nr);
            m2_avg=mean(m2_fit);
            m2_std=std(m2_fit);
            
            % Save results in the output data structure.
            S(j,k).Pfull(i)=log10(p_full);
            S(j,k).Pmlrg(i)=log10(p_mlrg);
            S(j,k).Mmax(i)=m2_avg;
            S(j,k).Mmax_err(i)=m2_std;
        end
    end
end

% Save the data.
save('TestBench_temp.mat','S');


