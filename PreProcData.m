function D=PreProcData(Cases)
  % Function that will pre-process the catalogue and injection information.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Get the data table and define an output strucutre.
  table_CaseDataTable;
  D=struct('Case',[],'Lat',[],'Lon',[],'Dep',[],'T',[],'M',[],'Mc',[],'Mk',[],'t',[],'v',[],'V',[]);
  
  % Loop over all the user-input cases.
  for i=1:length(Cases)
      
      % Get the table entry for this case.
      Case=Cases{i};
      D(i).Case=Case;
      j=find(strcmpi(Case,{T.Name}));
      source=T(j).source;
      load(T(j).Filename);
      D(i).Mc=T(j).Mc;
      D(i).Mk=T(j).Mk;
      
      % Process data to be the same format, depending on source.
      if(strcmpi(source,'Schultz22')) % From IS-Bath.
          % Catalogue information.
          D(i).Lat=S.cat.lat;
          D(i).Lon=S.cat.lon;
          D(i).Dep=S.cat.dep;
          D(i).T=S.cat.time;
          D(i).M=S.cat.mag;
          
          % Injection information.
          D(i).t=S.inj.time; % days.
          D(i).v=S.inj.rate; % m³/min.
          
          % Compute cumulative volume (m³).
          f=60*24; % minutes/day
          dt=D(i).t(2)-D(i).t(1);
          D(i).V=cumsum(D(i).v)*dt*f;
          
      elseif(strcmpi(source,'Study')) % From one off study.
          Dep
          
      elseif(strcmpi(source,'Verdon24')) % From Verdon23.
          % Catalogue information.
          D(i).Lat=[];
          D(i).Lon=[];
          D(i).Dep=[];
          D(i).T=R.T;
          D(i).M=R.M;
          
          % Injection information.
          D(i).t=R.IT; % days.
          D(i).v=R.IV; % m³/min.
          
          % Compute cumulative volume (m³).
          f=60*24; % minutes/day
          dt=D(i).t(2)-D(i).t(1);
          D(i).V=cumsum(D(i).v);

      end
  end
end


