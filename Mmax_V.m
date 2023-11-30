function [Mmax]=Mmax_V(V,params,type_flag)
  % Simple function that will compute maximum magnitudes for various models.
  %
  % V         - Cumulative volume injected (units of m^3).
  % params    - Vector of parameters.
  % type_flag - flag for the Mmax model to use.
  %
  % References:
  % 
  % Hanks, T. C., and H. Kanamori (1979), A moment magnitude scale, Journal of Geophysical Research, 84, 2348?2350, doi:10.1029/JB084iB05p02348.
  % Galis, M., Ampuero, J. P., Mai, P. M., & Cappa, F. (2017). Induced seismicity provides insight into why earthquake ruptures stop. Science advances, 3(12), eaap7528, doi: 10.1126/sciadv.aap7528.
  % McGarr, A. (1976). Seismic moments and volume changes. Journal of Geophysical Research, 81(8), 1487-1494, doi: 10.1029/JB081i008p01487.
  % McGarr, A. (2014). Maximum magnitude earthquakes induced by fluid injection. Journal of Geophysical Research: solid earth, 119(2), 1008-1019, doi: 10.1002/2013JB010597.
  % Van der Elst, N. J., Page, M. T., Weiser, D. A., Goebel, T. H., & Hosseini, S. M. (2016). Induced earthquake magnitudes are as large as (statistically) expected. Journal of Geophysical Research: Solid Earth, 121(6), 4575-4590, doi: 10.1002/2016JB012818.
  %
  % Written by Ryan Schultz.
  % 
  
  % The McGarr relationship [McGarr, 1976; 2014].
  if(strcmpi(type_flag,'McGarr'))
      
      % Fit parameter of apparent Shear Modulus (GPa).
      G=params(1);
      % Compute McGarr upperbound seismic moment (Eqn 13).
      Mo=G*V*1e9; % Seismic moment (N m).
      % and convert to moment magnitudes [Hanks & Kanamori, 1979].
      Mmax=(log10(Mo)-9.1)/1.5; % Seismic moment relationship using N m.

  % The van der Elst relationship [van der Elst et al., 2016].
  elseif(strcmpi(type_flag,'van der Elst'))
      
      % Fit parameter of Seismogenic Index and b-value.
      sigma=params(1);
      b=params(2);
      % Compute van der Elst Mmax (Eqns 4 & 5).
      Mmax=(log10(V)+sigma)/b;

  % The Galis relationship [Galis et al., 2017].
  elseif(strcmpi(type_flag,'Galis'))
      
      % Fit parameter of gamma (N m^(-7/2)).
      gamma=params(1);
      % Compute Galis maximum arrested rupture seismic moment (Eqn 3).
      Mo=gamma*(V.^1.5); % Seismic moment (N m).
      % and convert to moment magnitudes [Hanks & Kanamori, 1979].
      Mmax=(log10(Mo)-9.1)/1.5; % Seismic moment relationship using N m.

  end
  
end