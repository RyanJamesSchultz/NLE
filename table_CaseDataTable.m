% Simple script to tabulate all of the case data I'd like to examine in a structure.

% Define the strucutre table.
T=struct('Name',[], 'Filename',[], 'Mc',[], 'Mk',[], 'source',[]);
i=1;

% Basel.
T(i).Name='Basel';
T(i).Filename='data/BASEL_cat.mat';
T(i).Mc=1.5;
T(i).Mk=1.5;
T(i).source='Schultz22';
i=i+1;

% Cooper Basin, Habanero 1.
T(i).Name='CB1';
T(i).Filename='data/CB1_cat.mat';
T(i).Mc=0.0;
T(i).Mk=1.6;
T(i).source='Schultz22';
i=i+1;

% Cooper Basin, Habanero 4.
T(i).Name='CB4';
T(i).Filename='data/CB4_cat.mat';
T(i).Mc=2.0;
T(i).Mk=2.0;
T(i).source='Schultz22';
i=i+1;

% Paralana
T(i).Name='Paralana';
T(i).Filename='data/PARALANA_cat.mat';
T(i).Mc=0.5;
T(i).Mk=1.5;
T(i).source='Schultz22';
i=i+1;

% Pohang.
T(i).Name='Pohang';
T(i).Filename='data/POHANG_cat.mat';
T(i).Mc=0.5;
T(i).Mk=0.5;
T(i).source='Schultz22';
i=i+1;

% Duvernay ESB-10.
T(i).Name='ESB10';
T(i).Filename='data/ESB10_cat.mat';
T(i).Mc=0.8;
T(i).Mk=0.8;
T(i).source='Schultz22';
i=i+1;

% UK-PNR2.
T(i).Name='PNR2';
T(i).Filename='data/PNR2_cat.mat';
T(i).Mc=-1.0;
T(i).Mk=-0.5;
T(i).source='Schultz22';
i=i+1;

% Peace River.
T(i).Name='PeaceRiver';
T(i).Filename='data/PEACERIVER_cat.mat';
T(i).Mc=2.0;
T(i).Mk=2.0;
T(i).source='Schultz22';
i=i+1;

% Paradox Valley.
T(i).Name='ParadoxValley';
T(i).Filename='data/PARADOX_cat.mat';
T(i).Mc=1.5;
T(i).Mk=2.0;
T(i).source='Schultz22';
i=i+1;

% Guy-Greenbrier.
T(i).Name='GGB';
T(i).Filename='data/GUY_cat.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).source='Verdon24';
i=i+1;

% Prague.
T(i).Name='Prague';
T(i).Filename='data/PRAGUE_cat.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).source='Verdon24';
i=i+1;

% Pawnee.
T(i).Name='Pawnee';
T(i).Filename='data/PAWNEE_cat.mat';
T(i).Mc=2.5;
T(i).Mk=2.5;
T(i).source='Verdon24';
i=i+1;

% Youngstown.
T(i).Name='Youngstown';
T(i).Filename='data/YOUNGSTOWN_cat.mat';
T(i).Mc=1.4;
T(i).Mk=1.4;
T(i).source='Verdon24';
i=i+1;









% SSFS2000.
T(i).Name='SSFS2000';
T(i).Filename='data/SSFS2000.mat';
T(i).Mc=1.5;
T(i).Mk=1.5;
T(i).source='Schultz22';
i=i+1;

