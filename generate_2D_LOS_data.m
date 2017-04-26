function [y_POS,z_POS,resolution_meters,data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,data,r_obs,theta_obs,lambda_obs,FOV,resolution,u,G_factor)
%GENERATE_2D_LOS_DATA
%
% function [x_POS,y_POS,resolution_meters,data_2D_LOS] =
% generate_2D_LOS_data(y_data,z_data,z_data,data,r_obs,theta_obs,lambda_obs,FOV,resolution,u,G_factor)
% ...
%   *_data the grid CUBE...
%
%  data_2D_LOS is main output . it is 2D coro Im of True structure. ie the
%  synthetic images like EEGL.

%   Antti Pulkkinen, March 2017.

% Solar radius.
Rs = 695700e3; % m.
% 1 AU
AU = 149598000e3; % m.

% Resolution in radians.
resolution_rad = resolution*(pi/648000);
% Resolution in distance at 1 AU. https://en.wikipedia.org/wiki/Angular_diameter.
resolution_meters = tan(resolution_rad/2)*2*AU;

% Generate the LOS grid in the plane of sky (yz-plane).
[y_POS_tmp,z_POS_tmp] = meshgrid(-FOV(2)*Rs:resolution_meters:FOV(2)*Rs,-FOV(2)*Rs:resolution_meters:FOV(2)*Rs);

% Remove LOS points outside the FOV. NPS:creates the circle around the square image
radius_POS = sqrt(y_POS_tmp.^2 + z_POS_tmp.^2);
kk = find(radius_POS >= FOV(1)*Rs & radius_POS <= FOV(2)*Rs); y_POS = y_POS_tmp(kk); z_POS = z_POS_tmp(kk); x_POS = zeros(size(y_POS));

% Rotation matrices for rotation about z- and y-axis (theta and lambda, respectively).
theta_rad = theta_obs*pi/180; lambda_rad = lambda_obs*pi/180;
R_z = [cos(theta_rad) -sin(theta_rad) 0 ; sin(theta_rad) cos(theta_rad) 0 ; 0 0 1];
R_y = [cos(-lambda_rad) 0 sin(-lambda_rad) ; 0 1 0 ; -sin(-lambda_rad) 0 cos(-lambda_rad)];

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

% Grid definitions for the ray tracing algorithm.
grid3D.nx = length(x_data(1,:,1));
grid3D.ny = length(y_data(:,1,1));
grid3D.nz = length(z_data(1,1,:));
grid3D.minBound = [min(x_data(1,:,1)), min(y_data(:,1,1)), min(z_data(1,1,:))]';
grid3D.maxBound = [max(x_data(1,:,1)), max(y_data(:,1,1)), max(z_data(1,1,:))]';

% Initialize the data_2D_LOS. Default is zero emissions detected in the LOS element.
data_2D_LOS = zeros(1,length(x_POS_rot));

% Loop over the LOS.
for ii = 1:length(x_POS_rot),
    
    % LOS unit vector.
    r_LOS = -[(x_obs_rot - x_POS_rot(ii)) ; (y_obs_rot - y_POS_rot(ii)) ; (z_obs_rot - z_POS_rot(ii))];
    
    % Determine the grid indices the LOS ray passes through. THIS IS
    % WHERE WE CAN PLUG IN ALSO OTHER RAY TRACING ALGORITHMS.
    grid_indices = amanatidesWooAlgorithm_AP(r_obs_rot, r_LOS, grid3D); %pause(0.5); close;
    
    % Continue only if the LOS pierced the data cube through more than one grid point.
    if length(grid_indices) > 2,
        
        % Extract data values along the LOS. Values outside the data cube domain are set to 0.
        data_along_LOS = data(grid_indices);
        
        % Thomson scattering parameters along the LOS.
        [G_T,G_P,G_R,G_tot] = G_Thomson(x_obs_rot,y_obs_rot,z_obs_rot,x_data(grid_indices),y_data(grid_indices),z_data(grid_indices),u);
        
        % Segment lengths of the grid points along the LOS ray. ASSUME THE GRID INDICES ARE GIVEN IN ORDER FROM START TO END OF THE RAY.
        % YOU MAY WANT TO CHECK THIS MORE!
        cube_pierce_length = sqrt( (x_data(grid_indices(end)) - x_data(grid_indices(1))).^2 + (y_data(grid_indices(end)) - y_data(grid_indices(1))).^2 + (z_data(grid_indices(end)) - z_data(grid_indices(1))).^2 );
        dr_LOS_ray = cube_pierce_length/(length(grid_indices) - 1);
        
        % LOS integral.
        eval(sprintf('data_2D_LOS(ii) = sum(dr_LOS_ray*%s.*data_along_LOS);',G_factor));
        
    end;
    
end;
