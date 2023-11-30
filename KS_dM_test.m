function [KSp]=KS_dM_test(M,m1,resamp_flag)
  % Function that will use the Kolmogorov–Smirnov test to check for the 
  % influence of a maximum magnitude in a user-input catalogue by observing 
  % the sequencing of the largest magnitude events.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Make M a row vector and get its length.
  M=M(:)';
  Nm=length(M);
  %Nl=harmonic(Nm); % Approximated by log(Nm)+0.577215.
  
  % Get the sequence of differences in largest event magnitudes.
  [Mlrg,~]=OrderStatistic(M,Nm,'unique'); Mlrg=[m1,Mlrg];
  dMlrg=diff(Mlrg);
  
  % Do the M-permutation resampling of dMlrg.
  if(strcmpi(resamp_flag,'resample'))
      c=100;
      F=c; R=c/10;
      n=round(fminbnd(@(n) abs(F*harmonic(n)/n-R), 1,Nm));
      for j=1:F
          I=randperm(Nm); i=1;
          while(i<Nm)
              Ii=I(i:min(i+n-1,Nm));
              Mi=M(Ii);
              dmi=OrderStatistic(Mi,length(Mi),'unique'); dmi=[m1,dmi]; dmi=diff(dmi);
              dMlrg=[dMlrg,dmi];
              i=i+n;
          end
      end
  end
  
  % Compute the KS-test p-value.
  [~,KSp]=kstest2(M-m1,dMlrg);
  
end