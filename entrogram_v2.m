%
% CALCULATION OF THE ENTROGRAM (v1.0 - 24.02.2022)
%
function [H0,HL] = entrogram_v2(xyz_ind_gslib,Lagsin,bandwidth_X,bandwidth_Y,bandwidth_Z,...
    alpha,beta,NI,NJ,NK,Ncat,NMC)

% Calculate the global H--------------------------------------------------------
for n = 1:Ncat
    ind_all(:,n) = eq(xyz_ind_gslib,n);
end
pnow_all=sum(ind_all)./numel(xyz_ind_gslib);
lnpnow_all=log(pnow_all);
plnp_all=pnow_all.*lnpnow_all;
H0=-nansum(plnp_all);
%keyboard()
clear ind_all

% Start entrograms calculation

% Generate sampling locations using Sobol quasi random numbers
randpath = sobol_dataset(3,NMC,randi(1))';
count_lag=0;
for lag = Lagsin %this is the loop on lags
    %%
    if lag<max([bandwidth_X,bandwidth_Y,bandwidth_Z])
        count_lag=count_lag+1;
        textdisp= ['lag num = ',num2str(count_lag)];
        disp(textdisp)
        count_nmc = 0;
        %check lag vs bandwith
        if lag> bandwidth_X
            lag_x = bandwidth_X;
        else
            lag_x = lag;
        end
        if lag> bandwidth_Y
            lag_y = bandwidth_Y;
        else
            lag_y = lag;
        end
        if lag> bandwidth_Z
            lag_z = bandwidth_Z;
        else
            lag_z = lag;
        end
        for nmc=1:NMC %this is the loop on repetitions to compute for each lag
            % pick a node
            flag_dim = 0;
            ix = ceil(randpath(nmc,1)*NJ); %J index of the node
            iy = ceil(randpath(nmc,2)*NI); %I index of the node
            iz = ceil(randpath(nmc,3)*NK); %K index of the node
            % Define sector to search and check dimensions
            start_x = ix-lag_x;
            end_x = ix+lag_x;
            start_y = iy-lag_y;
            end_y = iy+lag_y;
            start_z = iz-lag_z;
            end_z = iz+lag_z;
            %Define all the points within the search sector
            X=start_x:end_x;
            Y=start_y:end_y;
            Z=start_z:end_z;
       
            [x,y,z] = meshgrid(X,Y,Z); %coordinates of all the pts within sector
            
            %Rotate around the z axis (alpha = strike = yaw, anticlockwise) 
            % and around the y axis   (beta = dip = pitch)
            xrotzy = ix + round((x-ix)*cos(deg2rad(alpha))*cos(deg2rad(beta))...
                -(y-iy)*sin(deg2rad(alpha))+(z-iz)*cos(deg2rad(alpha))*sin(deg2rad(beta)));
            yrotzy = iy + round((x-ix)*sin(deg2rad(alpha))*cos(deg2rad(beta))+...
                (y-iy)*cos(deg2rad(alpha))+(z-iz)*sin(deg2rad(alpha))*sin(deg2rad(beta)));
            zrotzy = iz + round(-(x-ix)*sin(deg2rad(beta))+(z-iz)*cos(deg2rad(beta)));
            
            
            % Check on the dimesions: if the sector boundaries fall outside
            % the domain, the sector is flagged and it is not considered 
            % for entropy calculation
            if ( min(min(min(xrotzy)))<1 || max(max(max(xrotzy)))>NJ )
                flag_dim = 1;
            end
            if ( min(min(min(yrotzy)))<1 || max(max(max(yrotzy)))>NI )
                flag_dim = 1;
            end
            
            if NK>1
                if ( min(min(min(zrotzy)))<1 || max(max(max(zrotzy)))>NK )
                    flag_dim = 1;
                end
            else
                zrotzy=1;
            end
           
            %calculate entrogram only if point is accepted
            if flag_dim == 0
                count_nmc =count_nmc+1; %counter
                %Find the GSLIB indexes
                ind_matx = xrotzy +(yrotzy-1)*NJ + (zrotzy-1)*NJ*NI;
                %reshape them into array
                sector_ind = reshape(ind_matx,1,length(X)*length(Y)*length(Z));
                % Calculate local entropy
                class_sect =xyz_ind_gslib(sector_ind)'; %These is the array of classes within the search sector         for n = 1:Ncat
                ind_sect = zeros(numel(class_sect),Ncat);
                for n = 1:Ncat
                    ind_sect(:,n)=eq(class_sect,n);
                end
                pnow_tmp=squeeze(sum(ind_sect))/((numel(X)*numel(Y)*numel(Z)));
                lnpnow=log(pnow_tmp);
                plnp_local=(pnow_tmp.*lnpnow)';
                HL_temp(count_nmc) =-nansum(plnp_local);
            end
        end
        %Calculate the unit lag distance for the entrogram
        if count_nmc>0
            %Magnitude of vector
            %dist(count_lag,1) = sqrt((lag_x*dx)^2+(lag_y*dy)^2+(lag_z*dz)^2);
            HL(count_lag,1)=mean(HL_temp); %entrogram
            clear HL_temp;
        end
    else
        break
    end
end


