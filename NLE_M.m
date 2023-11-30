function [Mq]=NLE_M(m1,m2,b,q)
  % Function to determine the magnitude of the next largest event (NLE) 
  % in a sequence, given the input catalogue and multiple model estimates 
  % of Mmax.
  % Code is vectorized.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Useful constant.
  cn=1-10.^(-b.*(m2-m1));
  
  % Flag for output type.
  if(any(q==-1))

      % Use the analytical estimate of the mean NLE magnitude.
      cd=b*log(10);
      if(any(isinf(m2)))
          Mq=(m1.*cd+1)./cd;
      else
          Mq=((m1.*cd+1)-(1-cn).*(m2.*cd+1))./(cd.*cn);
      end
      
  else
      
      % Use an inverted GR-MFD CDF to map the quantile q into magnitudes.
      Mq=m1-log10(1-q.*cn)./b;
  end
  
end

