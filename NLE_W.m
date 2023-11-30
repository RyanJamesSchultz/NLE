function [W]=NLE_W(M,dMlrg,Mlrg,m1,m2,b)
  % Function to determine the ensemble weights for next largest event 
  % (NLE), given a sequence dMlrg (with Mlrg) and multiple model estimates 
  % of Mmax.
  % 
  % References:
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Preallocate.
  LLR=zeros(size(m2));
  AIC=LLR; BIC=LLR;
  
  % Loop over all of the Mmax models.
  for j=1:length(m2)
      
      % Get the AIC/BIC scores for both the full and dM catalogues.
      [LLn,AICn,BICn]=GR_MFD_LL(    M, m1, m2(j),     0,b); AICn=AICn-2*LLn; BICn=BICn-2*LLn;
      [LLd,AICd,BICd]=GR_MFD_LL(dMlrg,  0, m2(j)-Mlrg,0,b); AICd=AICd-2*LLd; BICd=BICd-2*LLd;
      
      % Error handling.
      if(AICn==Inf)
          AICn=BICn; % This is only fine if K is constant across models.
      end
      if(AICd==Inf)
          AICd=BICd; % This is only fine if K is constant across models.
      end
      
      % Combine scores.
      LLR(j)=(LLd-LLn)*2;
      AIC(j)=(AICn-AICd);
      BIC(j)=(BICn-BICd);
      
  end
  
  % Get the ensemble weights [Wagenmakers & Farrell, 2004].
  dAIC=AIC-min(AIC); Wa=exp(-dAIC/2); Wa=Wa/sum(Wa);
  dBIC=BIC-min(BIC); Wb=exp(-dBIC/2); Wb=Wb/sum(Wb);
  W=Wa+Wb; W=W/sum(W);
  
end

