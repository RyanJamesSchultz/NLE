function [x,w]=meanNLL(X,nLL,p)
  % Function that computes the weighted mean of a sample X, given the input
  % negative log-likelihood scores (e.g., AIC or BIC).
  % 
  % Written by Ryan Schultz.
  % 
  
  % Unwrap X & nLL into a column vector.
  X=X(:);
  nLL=nLL(:);
  
  % Get the weights.
  dLL=nLL-min(nLL);
  w=exp(-dLL/2);
  w=w/sum(w);
  
  % Discard the bottom p-th percentile of weights.
  w(w<prctile(w,p))=0;
  w=w/sum(w);
  
  % Get the weighted mean.
  x=sum(w.*X);
  
end