function [S]=EnsembleW(M,m1,S_m2,b,smooth_flag)
  % Function to determine the sequence of ensemble weights for next largest 
  % events (NLEs), given the input catalogue and multiple model estimates 
  % of Mmax.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Make M a row vector and get some important lengths.
  M=M(:)';
  Nc=length(M);
  Nm=length(S_m2);
  
  % Get the sequence of largest events (and indicies to when they occur).
  Mlrg=OrderStatistic(M,Nc,'none'); Mlrg=[m1,Mlrg];
  [Mlrg,I]=unique(Mlrg);
  dMlrg=diff(Mlrg);
  Mlrg(end)=[];
  
  % Predefine/preallocate the output weight structure.
  S=struct('Mmax',[],'W',[]);
  for j=1:Nm
      S(j).W=ones([1 Nc+1])/Nm;
      if(length(S_m2(j).Mmax)==(Nc+1))
          S(j).Mmax=S_m2(j).Mmax;
      else
          S(j).Mmax=S_m2(j).Mmax*ones([1 Nc+1]);
      end
  end
  wp=ones([1 Nm])/Nm;
  
  % Loop over all of the unique Mlrg entries.
  for k=2:length(I)
      
      % Get the ensemble weights for all of the Mlrg values up to now.
      i=I(k);
      Mmax=arrayfun(@(S) S.Mmax(i),S);
      w=NLE_W(M(1:i-1),dMlrg(1:k-1),Mlrg(1:k-1),m1,Mmax,b);
      
      % Optionally, smooth out the sequence of weights.
      if(strcmpi(smooth_flag,'mean'))
          w=w+wp;
          w=w/sum(w);
      elseif(strcmpi(smooth_flag,'geomean'))
          w=sqrt(w.*wp);
          w=w/sum(w);
      elseif(strcmpi(smooth_flag,'both'))
          w1=w+wp;        w1=w1/sum(w1);
          w2=sqrt(w.*wp); w2=w2/sum(w2);
          w=w1+w2;
          w=w/sum(w);
      end
      wp=w;
      
      % Stuff the results into the output data structure.
      for j=1:Nm
          S(j).W(i:end)=w(j);
      end
  end
  
end

