function [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = map_LOS_2_G_data(x_data,y_data,z_data,y_POS,z_POS,r_obs,theta_obs,lambda_obs,u)
%MAP_LOS_2_G_DATA
%
% [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = generate_LOS_2_G_data(x_data,y_data,z_data,y_POS,z_POS,r_obs,theta_obs,lambda_obs,u) ...
%
%   *_data the grid CUBE...
%
%  *_POS - plane-of-sky coordinates in the synthetic coronagraph image.
%
% Note: it is assumed that the grid is a cube. Other grid shapes will lead
% to erroneous LOS tracking within the volume.

%   Antti Pulkkinen, March 2017.

% Rotation matrices for rotation about z- and y-axis (theta and lambda, respectively).
theta_rad = theta_obs*pi/180; lambda_rad = lambda_obs*pi/180;
R_z = [cos(theta_rad) -sin(theta_rad) 0 ; sin(theta_rad) cos(theta_rad) 0 ; 0 0 1];
R_y = [cos(-lambda_rad) 0 sin(-lambda_rad) ; 0 1 0 ; -sin(-lambda_rad) 0 cos(-lambda_rad)];

x_POS = zeros(size(y_POS));

% Rotate the plane of sky to correspond to the location of the observer.
x_POS_rot = zeros(size(x_POS)); y_POS_rot = x_POS_rot; z_POS_rot = x_POS_rot;
for ii = 1:length(x_POS),
    
    r_POS_rot_zy = R_z*R_y*[x_POS(ii) ; y_POS(ii) ; z_POS(ii)];
    
    x_POS_rot(ii) = r_POS_rot_zy(1);
    y_POS_rot(ii) = r_POS_rot_zy(2);
    z_POS_rot(ii) = r_POS_rot_zy(3);
    
end;

% Rotate the observer coordinates.
x_obs = r_obs; y_obs = 0; z_obs = 0;
r_obs_rot = R_z*R_y*[x_obs ; y_obs ; z_obs];
x_obs_rot = r_obs_rot(1); y_obs_rot = r_obs_rot(2); z_obs_rot = r_obs_rot(3);

% Loop over the LOS vectors.
% NPS: This FOR LOOP, is repeat subroutine from generate_2D_LOS_data
% NPS: I am working on making this code work as a single line replacement

% OUT_Indices = NPS_eg_bresen (Ist_loc, Ien_loc, IN_extra);

% Grid definitions for the ray tracing algorithm.
grid3D.nx = length(x_data(1,:,1));
grid3D.ny = length(y_data(:,1,1));
grid3D.nz = length(z_data(1,1,:));
grid3D.minBound = [min(x_data(1,:,1)), min(y_data(:,1,1)), min(z_data(1,1,:))]';
grid3D.maxBound = [max(x_data(1,:,1)), max(y_data(:,1,1)), max(z_data(1,1,:))]';

% Maximum possible number of grid point in a given LOS. This will be
% used to initialize the grid index and other matrices.
max_no_grid_points = ceil(sqrt( grid3D.nx.^2 + grid3D.ny.^2 + grid3D.nz.^2 ));

grid_indices = NaN*zeros(length(y_POS),max_no_grid_points); cube_pierce_length = NaN*zeros(length(y_POS),1);
G_T_LOS = grid_indices; G_P_LOS = grid_indices; G_R_LOS = grid_indices; G_tot_LOS = grid_indices;

% Loop over LOS (i.e. plane-of-sky points)
for ii = 1:length(x_POS_rot),
    
    % LOS vector.
    r_LOS = -[(x_obs_rot - x_POS_rot(ii)) ; (y_obs_rot - y_POS_rot(ii)) ; (z_obs_rot - z_POS_rot(ii))];
    
    % Determine the grid indices that the LOS ray passes through. THIS IS
    % WHERE WE CAN PLUG IN ALSO OTHER RAY TRACING ALGORITHMS.
    grid_indices_LOS = amanatidesWooAlgorithm_AP(r_obs_rot, r_LOS, grid3D);
    
    % Check if the LOS went through the data domain via more than one grid point. If not, return NaNs.
    if length(grid_indices_LOS) > 1,
        
        % Segment lengths of the grid points along the LOS ray. ASSUME THE GRID INDICES ARE GIVEN IN ORDER FROM START TO END OF THE RAY.
        cube_pierce_length(ii) = sqrt( (x_data(grid_indices_LOS(end)) - x_data(grid_indices_LOS(1))).^2 + (y_data(grid_indices_LOS(end)) - y_data(grid_indices_LOS(1))).^2 + (z_data(grid_indices_LOS(end)) - z_data(grid_indices_LOS(1))).^2 );
                        
        % Record to the master indices holder.
        grid_indices(ii,1:length(grid_indices_LOS)) = grid_indices_LOS;
                
%         % Coordinates at which G-factors are evaluated.
%         x_Gfactor = x_data(grid_indices(ii,isfinite(grid_indices(ii,:))));
%         y_Gfactor = y_data(grid_indices(ii,isfinite(grid_indices(ii,:))));
%         z_Gfactor = z_data(grid_indices(ii,isfinite(grid_indices(ii,:))));

        % Coordinates at which G-factors are evaluated.
        x_Gfactor = x_data(grid_indices_LOS);
        y_Gfactor = y_data(grid_indices_LOS);
        z_Gfactor = z_data(grid_indices_LOS);
        
        % Thomson scattering parameters along the LOS at the grid points.
        [G_T_tmp,G_P_tmp,G_R_tmp,G_tot_tmp] = G_Thomson(x_obs_rot,y_obs_rot,z_obs_rot,x_Gfactor,y_Gfactor,z_Gfactor,u);
        
        G_T_LOS(ii,1:length(G_T_tmp)) = G_T_tmp;
        G_P_LOS(ii,1:length(G_P_tmp)) = G_P_tmp;
        G_R_LOS(ii,1:length(G_R_tmp)) = G_R_tmp;
        G_tot_LOS(ii,1:length(G_tot_tmp)) = G_tot_tmp;
        
    end;
    
end;
