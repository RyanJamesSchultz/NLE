function [Mk,I]=OrderStatistic(M,k,truncation_flag)
  % Function that computes the (time-dependent) order statisitc of some 
  % sample.  This function assumes that the input sample is already  
  % chronologically sorted.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the sample size.
  N=length(M);
  I=[];
  
  % If-block to choose the algorithm that will find the order statistic.
  if(k==1) % Finding the smallest value.
      Mk=cummin(M);
      
  elseif(k==N) % Finding the largest value.
      Mk=cummax(M);
      
  else % Anything else in between.
      
      % Start looping from largest value and reduce until the k-th order statistic.
      for i=N:-1:k
          % Select the current largest values in sample.
          Mi=cummax(M);
          
          % and remove them from the sample.
          [~,Il]=unique(Mi);
          M(Il)=NaN; % This will produce NaN values when there isn't a k-th order statistic yet.
      end
      Mk=Mi;
      
  end
  
  % Only keep the (first) unique record values, if flagged to.
  if(strcmpi(truncation_flag,'unique'))
      [Mk,I]=unique(Mk);
  end
  
end
