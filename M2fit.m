function [m2,LLR]=M2fit(M,m1,b,Nr)
  % Function that will MLE fit the Mmax value, using magnitude differences 
  % in the sequence of the largest events.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Make M a row vector and get the catalogue length.
  M=M(:)';
  Nm=length(M);
  
  % Get the sequence of differences in largest event magnitudes.
  [Mlrg,~]=OrderStatistic(M,Nm,'unique'); Mlrg=[m1,Mlrg];
  dMlrg=diff(Mlrg);
  Mlrg(end)=[];
  
  % Anonymous function that computes the (unscaled) log-likelihood ratio.
  LR = @(m) +2*GR_MFD_LL(dMlrg,0,m-Mlrg,0,b)-2*GR_MFD_LL(M,m1,m,0,b);
  
  % Search for the optimal m2 parameter using a Maximum Likelihood Estimator (MLE).
  options = optimoptions(@fmincon,'StepTolerance',1e-4,'ConstraintTolerance',1e-3,'Display','none');
  m2 = fmincon(LR,max(M)+0.2,[],[],[],[],max(M),10,[],options);
  
  % Compute the (scaled) log-likelihood ratio.
  LLR=+2*GR_MFD_LL(dMlrg,0,m2-Mlrg,0,b)-2*GR_MFD_LL(M,m1,m2,0,b)-2*GR_MFD_LL(dMlrg,0,Inf,0,b)+2*GR_MFD_LL(M,m1,Inf,0,b);
  
  % Do the M-permutation resampling of dMlrg.
  for j=2:Nr+1
      
      % Get a new catalogue by reshuffling.
      Ii=randperm(Nm); Mi=M(Ii);
      
      % Get a new Mlrg & dMlrg sequence.
      Mlrg_i=OrderStatistic(Mi,length(Mi),'unique'); Mlrg_i=[m1,Mlrg_i]; dMlrg_i=diff(Mlrg_i); Mlrg_i(end)=[];
      
      % Find the optimal m2 value (and the corresponding log-likelihood ratio).
      LR = @(m)  +2*GR_MFD_LL(dMlrg_i,0,m-Mlrg_i,0,b)-2*GR_MFD_LL(Mi,m1,m,0,b);
      m2(j)=fmincon(LR,max(M)+0.2,[],[],[],[],max(M),10,[],options);
      LLR(j)=+2*GR_MFD_LL(dMlrg_i,0,m2(j)-Mlrg_i,0,b)-2*GR_MFD_LL(Mi,m1,m2(j),0,b)-2*GR_MFD_LL(dMlrg_i,0,Inf,0,b)+2*GR_MFD_LL(Mi,m1,Inf,0,b);
  end
  
end 
