function [LL,AICc,BIC]=GR_MFD_LL(M,m1,m2,a,b)
  % Function that computes the log-likelihood (and similar scores) for the 
  % Gutenberg-Richter magnitude-frequency distributions (GR-MFD), given a 
  % sample catalogue.
  % 
  % References:
  % Aki (1965). Maximum likelihood estimate of b in the formula log (N) = a − bM and its confidence limits, Bull. Earthq. Res. Inst. Tokyo Univ., 43, 237-239.
  % Marzocchi & Sandri (2003). A review and new insights on the estimation of the b-valueand its uncertainty. Annals of geophysics.
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the upper magnitude bounds and length.
  n=length(M);
  K=3+1; % b,m1,m2.
  
  % Get the GR-MFD log-likelihood [Aki, 1965; Marzocchi & Sandri, 2003].
  [PDF,~,~]=GR_MFD(M, m1,m2, a,b, 'norm');
  LL=sum(log(PDF));
  
  % Compute the AIC & BIC statistics [Wagenmakers & Farrell, 2004].
  %AIC=2*K-2*LL;
  AICc=2*K+(2*K*(K+1)/(n-K-1))-2*LL;
  BIC=K*log(n)-2*LL;
  
end