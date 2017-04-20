function IndexGrd = voxel2ind (lenx,leny,vox_xyz)

% lenx      = total voxels along x direction
% leny      = total voxels along y direction
% lenz      = total voxels along z direction
% vox_xyz   = [ix,iy,iz] - grid index of each axis. 
%                          Each axis must be a colm vector 
% 
% TEST SETUP:
% [xg,yg,zg] = meshgrid(0:0.5:12,0:1:24,0:2:48); 
% xgdpts = xg(1,:,1); ygdpts = yg(:,1,1)'; zgdpts = squeeze(zg(1,1,:))';
% entry_P = [5.5,3.1,0.2];
% lenx = length(xgdpts);
% leny = length(ygdpts);

% ix =  find( (xgdpts-entry_P(1,1))<=0,1,'last') ;
% iy =  find( (ygdpts-entry_P(1,2))<=0,1,'last') ;
% iz =  find( (zgdpts-entry_P(1,3))<=0,1,'last') ;
% vox_xyz = [ix,iy,iz];
% 
% IndexGrd = voxel2ind (lenx,leny,vox_xyz)
%
ix = vox_xyz(:,1);
iy = vox_xyz(:,2);
iz = vox_xyz(:,3);

IndexGrd = (lenx*leny)*(iz-1) + leny *(ix-1) + iy;
        
        
        
end