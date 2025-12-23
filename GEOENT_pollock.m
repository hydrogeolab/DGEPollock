clear all; close all; fclose all; clc;

bandwidth_X=50;         % max size (number of cells) of searching box along X (J) direction
Ncat=10;                % number of categories (nb)
thre=0;                 % resolution threshold
a=dir("*.jpg");

for j=1:numel(a)

    I= imread(a(j).name);
    IX(j)=size(I,1);
    IY(j)=size(I,2);

end
wHD=find(IX>thre);

a=a(wHD);

for j=1:numel(a)
    
    close all;
    c=strsplit(a(j).name,'_');

    yea(j)=str2num(c{1});

    LY(j)=str2num(c{2}); %cm
    LX(j)=str2num(c{3}); %cm

    Data=importdata(a(j).name);

    L=0.299*Data(:,:,1) + 0.587 * Data (:,:,2) + 0.114 * Data (:,:,3);
      
    NJ=size(Data,2);            % number of cells in X direction
    NI=size(Data,1);            % number of cells in Y direction
    NK=1;  

    field_3D=reshape(L,[NJ,NI,NK]); % Matlab 3D arrays: Y,X,Z
    field_3D=flipud(rot90(field_3D));

    dx = NJ/LX(j);            % size/spacingj of the cells in X direction
    dy = NI/LY(j);            % size/spacing of the cells in Y direction
    dz = 1.;            % size/spacing of the cells in Z direction

    input_gslib=vtk2gslib(field_3D,NI,NJ,NK);
    
    str1=a(j).name;

    close all

    NMC=200;                % number of repetitions to compute ensemble average of entrograms (min 2 - suggested 500)
    loglag=false;           % if active: lags are computed logarithmically
    lagsteplogspace=1000;   % if loglag=true: steps of the logarithmically spaced lags
    lagstep=1;              % if loglag=false: frequency of the linearly spaced lags

    alpha = [0. 90. 0.];    % angle of rotation of the searching box about the z axis ("strike")
    beta = [0. 0. 90.];     % angle of rotation of the searching box about the y axis ("dip")
    nHs_flag=1;             % if active: normalized entropic scale is reported

   
    bandwidth_Y=1;          % max size (number of cells) of searching box along Y (I) direction
    bandwidth_Z=1;          % max size (number of cells) of searching box along Z (K) direction (set 1 for 2D datasets)
    NUMDIR = 2;             % number of directions (max 3)
    
    GEOENT_core_v2
    
    H0_all_anis1(j)=H0_nd(1);
    H0_all_anis2(j)=H0_nd(2);
    nHs_all_anis1(j)=normHscale(1);
    nHs_all_anis2(j)=normHscale(2);
       
    filename=strcat(a(j).name,'_anis.fig');
    saveas(gcf,filename);
    filename=strcat(a(j).name,'_anis.png');
    exportgraphics(gcf,filename); % uncomment to save as PNG file.
    filename=strcat(a(j).name,'_anis.mat');
    save(filename,'-v7.3'); % uncomment to store results in .mat file.

    bandwidth_Y=bandwidth_X;      % max size (number of cells) of searching box along Y (I) direction
    NUMDIR = 1;          % number of directions (max 3)
    GEOENT_core_v2

    filename=strcat(a(j).name,'_iso.fig');
    saveas(gcf,filename);
    filename=strcat(a(j).name,'_iso.png');
    exportgraphics(gcf,filename); % uncomment to save as PNG file.
    filename=strcat(a(j).name,'_iso.mat');
    save(filename,'-v7.3'); % uncomment to store results in .mat file.

    close all

    H0_all_iso(j)=H0_nd(1);
    nHs_all_iso(j)=normHscale(1);
    HR_all_iso(j)=HL(1)/H0_nd(1);

end

