% A fast and simple voxel traversal algorithm through a 3D space partition (grid)
% proposed by J. Amanatides and A. Woo (1987).

cc;

% % Test Nro. 1
origin    = [15, 15, 15]';
direction = [-0.3, -0.5, -0.7]';

% % Test Nro. 2
%origin    = [-8.5, -4.5, -9.5]';
%direction = [0.5, 0.5, 0.7]';

[x,y,z] = meshgrid(-12:12,-15:15,-10:10);

% Grid: dimensions
grid3D.nx = length(x(1,:,1));
grid3D.ny = length(y(:,1,1));
grid3D.nz = length(z(1,1,:));
grid3D.minBound = [min(x(1,:,1)), min(y(:,1,1)), min(z(1,1,:))]';
grid3D.maxBound = [max(x(1,:,1)), max(y(:,1,1)), max(z(1,1,:))]';

% % Grid: dimensions
% grid3D.nx = 10;
% grid3D.ny = 15;
% grid3D.nz = 20;
% grid3D.minBound = [-10, -10, -20]';
% grid3D.maxBound = [ 10,  10,  20]';
 
verbose = 1;
voxel_indices = amanatidesWooAlgorithm_AP(origin, direction, grid3D, verbose);

figure;
plot3(x(:),y(:),z(:),'.'); grid on; hold on; plot3(x(voxel_indices),y(voxel_indices),z(voxel_indices),'o');
plot3(x(voxel_indices(1)),y(voxel_indices(1)),z(voxel_indices(1)),'.','markersize',50);
plot3(x(voxel_indices(end)),y(voxel_indices(end)),z(voxel_indices(end)),'.','markersize',50);
xlabel('x'); ylabel('y'); zlabel('z');
