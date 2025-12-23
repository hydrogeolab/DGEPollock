function gslibvector=import_raw(inputfilename,m_height, m_width, m_depth)

nz=m_height;
nx=m_width;
ny=m_depth;

VoxelMat= multibandread(inputfilename, [m_height,m_width,m_depth],'uint8',0,'bip','ieee-le');

% transform to GSLIB
gslibvector=zeros(1000^3,1);
for iz=1:nz
    for iy=1:ny
        for ix=1:nx
            loc=(iz-1)*nx*ny + (iy-1)*nx + ix;
            gslibvector(loc)= VoxelMat(iy,ix,iz);
        end
    end
end


% save('gslibfileIN.mat', 'gslibvector','-v7.3');
%& create a file to visualize it in SGEMS:
%dlmwrite('microCTgslib.dat',gslibvector(1:1000^3),'precision',0);
%fclose all;