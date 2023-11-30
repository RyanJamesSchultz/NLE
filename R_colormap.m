function [Z] = R_colormap(type_flag)
  % colormap()
  % 
  % Simple function to generate a color maps to display odds ratio information.
  %
  % Written by Ryan Schultz.
  %
  
  % Define the colormaps.
  if(strcmpi(type_flag,'Odds Ratio'))
      n=50;
      Z=[82      0     112;  ... % Dark magenta.
         255     0     225;  ... % Magenta.
         255     255   255;  ... % White.
         255     0     0  ;  ... % Red.
         107     13    14 ]; ... % Dark red.
     Z=flipud(Z);
  else
      n=50;
      Z=0;
  end
  
  % Make 16-bit RGB values between 0-1.
  Z=Z/255;
  
  % Interpolate to the desired number of samples.
  I=1:length(Z);
  Z=interp1(I,Z,min(I):1/n:max(I));
   
return;