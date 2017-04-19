function [y_POS,z_POS,resolution_meters,data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,data,r_obs,theta_obs,lambda_obs,FOV,resolution,dr,u,G_factor)
%GENERATE_2D_LOS_DATA
%
% function [x_POS,y_POS,resolution_meters,data_2D_LOS] =
% generate_2D_LOS_data(y_data,z_data,z_data,data,r_obs,theta_obs,lambda_obs,FOV,resolution,dr,u,G_factor)
% ...
%  data_2D_LOS is main output . it is 2D coro Im of True structure. ie the
%  syntheic images like EEGL.

%   Antti Pulkkinen, March 2017.

% Solar radius.
Rs = 695700e3; % m.
% 1 AU
AU = 149598000e3; % m. 

% Do not process LOS points outside the max. radius of the data domain.
LOS_limit_domain = max(sqrt( x_data(:).^2 + y_data(:).^2 + z_data(:).^2 )); % m.

% Resolution in radians.
resolution_rad = resolution*(pi/648000);
% Resolution in distance at 1 AU. https://en.wikipedia.org/wiki/Angular_diameter.
resolution_meters = tan(resolution_rad/2)*2*AU;

% Generate the LOS grid in the plane of sky (yz-plane).
[y_POS_tmp,z_POS_tmp] = meshgrid(-FOV(2)*Rs:resolution_meters:FOV(2)*Rs,-FOV(2)*Rs:resolution_meters:FOV(2)*Rs);

% Remove LOS points outside the FOV.
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

% Loop over the LOS vectors.
data_2D_LOS = NaN*x_POS_rot;
for ii = 1:length(x_POS_rot),
    
    % LOS unit vectors.
    r_LOS = [(x_obs_rot - x_POS_rot(ii)) ; (y_obs_rot - y_POS_rot(ii)) ; (z_obs_rot - z_POS_rot(ii))];
    e_LOS = r_LOS/sqrt(r_LOS.'*r_LOS);
    
    % Number of LOS steps required to march from the observer to the opposite side of the heliosphere.
    no_of_LOS_steps = length(-r_obs:dr:r_obs);
    
    r_LOS_steps = NaN*zeros(3,no_of_LOS_steps);
    radius_LOS_steps = NaN*zeros(1,no_of_LOS_steps);
    
    for jj = 1:no_of_LOS_steps;
        
        r_LOS_steps(:,jj) = r_obs_rot - jj*dr*e_LOS;
        radius_LOS_steps(:,jj) = sqrt(r_LOS_steps(:,jj).'*r_LOS_steps(:,jj));
    end;
    
    % For computational efficiency, include LOS points only within given LOS_limit_domain radius.
    kk_LOS_limit_domain = find(radius_LOS_steps < LOS_limit_domain);
    % Trim the LOS set.
    r_LOS_steps = r_LOS_steps(:,kk_LOS_limit_domain);
    
    % Interpolate data values along the LOS. Values outside the data cube domain are set to 0.
    data_along_LOS = interp3(x_data,y_data,z_data,data,r_LOS_steps(1,:),r_LOS_steps(2,:),r_LOS_steps(3,:),'nearest',0);
    
    % Thomson scattering parameters along the LOS.
    [G_T,G_P,G_R,G_tot] = G_Thomson(x_obs_rot,y_obs_rot,z_obs_rot,r_LOS_steps(1,:),r_LOS_steps(2,:),r_LOS_steps(3,:),u);
          
    % LOS integral.
    eval(sprintf('data_2D_LOS(ii) = sum(dr*%s.*data_along_LOS);',G_factor));
    
end;

% figure; plot3(x_POS_rot/Rs,y_POS_rot/Rs,z_POS_rot/Rs,'.'); grid on; xlabel('x'); ylabel('y'); zlabel('z'); axis equal;
% 
% figure; plot3(x_POS_rot/AU,y_POS_rot/AU,z_POS_rot/AU,'.'); grid on; xlabel('x'); hold on; ylabel('y'); zlabel('z'); axis equal;
% plot3(x_obs_rot/AU,y_obs_rot/AU,z_obs_rot/AU,'*'); xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
% plot3(r_LOS_steps(1,:)/AU,r_LOS_steps(2,:)/AU,r_LOS_steps(3,:)/AU,'color','k');
% 
% data_2D_LOS_grid = griddata(y_POS,z_POS,data_2D_LOS,y_POS_tmp,z_POS_tmp,'linear');
% [kk_i,kk_j] = find(radius_POS < FOV(1)*Rs); data_2D_LOS_grid(kk_i,kk_j) = NaN;
% figure; pcolor(y_POS_tmp/Rs,z_POS_tmp/Rs,data_2D_LOS_grid); colorbar;
% 
% pause;
% 
% 
% 
