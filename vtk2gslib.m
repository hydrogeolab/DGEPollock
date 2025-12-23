function xyz_ind_gslib=vtk2gslib(field_3D,NI,NJ,NK)
    
    xyz_ind_gslib=zeros(NI*NJ*NK,1);
    for iz=1:NK
        for iy=1:NI
            for ix=1:NJ    
                loc = (iz-1)*NI*NJ + (iy-1)*NJ + ix;
                xyz_ind_gslib(loc)=field_3D(iy,ix,iz);
            end
        end
    end
    