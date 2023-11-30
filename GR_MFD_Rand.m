function [M]=GR_MFD_Rand(m1,m2,a,b,Nr)
  % Function that randomly draws earthquake magnitudes from the 
  % Gutenberg-Richter magnitude-frequency distribution (GR-MFD) using
  % inverse transform sampling.
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % The PDF & CDF for a bounded GR-MFD.
  %cd=b*log(10);
  %x=10.^(-b*(m-m1));
  cn=1-10.^(-b.*(m2-m1));
  
  % Get the random CDF values.
  r=rand(Nr);
  
  % Use an inverted GR-MFD CDF to map the random CDF values into magnitudes.
  M=m1-log10(1-r.*cn)./b;
  
end
