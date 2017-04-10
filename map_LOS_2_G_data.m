function [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = map_LOS_2_G_data(x_data,y_data,z_data,y_POS,z_POS,r_obs,theta_obs,lambda_obs,dr,u)
%MAP_LOS_2_G_DATA
%
% [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = generate_LOS_2_G_data(x_data,y_data,z_data,y_POS,z_POS,r_obs,theta_obs,lambda_obs,dr,u) ...
%
% Note: dr should be chosen so that it is smaller than the grid spacing to
% ensure that all grid points are captured in the extraction of the G factors.
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

% Number of steps needed to go through the grid with given dr. This will be
% used to initialize the grid index and other matrices.
no_of_steps = round(max([( max(x_data(:)) - min(x_data(:))) (max(x_data(:)) - min(x_data(:))) (max(x_data(:)) - min(x_data(:)))])/dr);

grid_indices = NaN*zeros(length(y_POS),no_of_steps); cube_pierce_length = NaN*zeros(length(y_POS),1);
G_T_LOS = grid_indices; G_P_LOS = grid_indices; G_R_LOS = grid_indices; G_tot_LOS = grid_indices;

% Loop over the LOS vectors.
for ii = 1:length(x_POS_rot),
    
    % LOS unit vectors.
    r_LOS = [(x_obs_rot - x_POS_rot(ii)) ; (y_obs_rot - y_POS_rot(ii)) ; (z_obs_rot - z_POS_rot(ii))];
    e_LOS = r_LOS/sqrt(r_LOS.'*r_LOS);
    
    % Number of LOS steps required to ray trace from the observer to the opposite side of the heliosphere.
    no_of_LOS_steps = length(-r_obs:dr:r_obs);
    
    r_LOS_steps = NaN*zeros(3,no_of_LOS_steps);
    
    for jj = 1:no_of_LOS_steps;
        r_LOS_steps(:,jj) = r_obs_rot - jj*dr*e_LOS;    
    end;
    
    % Include LOS points only within the domain (assumed cube).
    kk_LOS_in_domain_xy = intersect(find(min(x_data(:)) <= r_LOS_steps(1,:) & r_LOS_steps(1,:) <= max(x_data(:))),find(min(y_data(:)) <= r_LOS_steps(2,:) & r_LOS_steps(2,:) <= max(y_data(:))));
    kk_LOS_in_domain = intersect(find(min(z_data(:)) <= r_LOS_steps(3,:) & r_LOS_steps(3,:) <= max(z_data(:))),kk_LOS_in_domain_xy);
    % Trim the LOS set.
    r_LOS_steps = r_LOS_steps(:,kk_LOS_in_domain);
    
    % Check if the LOS went through the data domain. If not, return NaNs.
    if ~isempty(kk_LOS_in_domain),
        
        closest_index = NaN*zeros(1,length(kk_LOS_in_domain));
        
        % Length of the LOS ray piercing the cube.
        cube_pierce_length(ii) = sqrt( (r_LOS_steps(1,1) - r_LOS_steps(1,end)).^2 + (r_LOS_steps(2,1) - r_LOS_steps(2,end)).^2 + (r_LOS_steps(3,1) - r_LOS_steps(3,end)).^2 );
        
        % Find the closest grid points along the LOS.
        for ss = 1:length(kk_LOS_in_domain),
            
            grid_distance = sqrt( (r_LOS_steps(1,ss) - x_data).^2 + (r_LOS_steps(2,ss) - y_data).^2 + (r_LOS_steps(3,ss) - z_data).^2 );
            [closest_grid_distance,closest_index(ss)] = min(grid_distance);
        end;
        
        % Make sure indices are not repeated and record to master indices holder.
        tmp_indices = unique(closest_index);
        grid_indices(ii,1:length(tmp_indices)) = tmp_indices;
        
        
%             %% TEST
%               % Solar radius.
%             Rs = 695700e3; % m.
% 
%  %          figure; plot3(x_POS_rot/Rs,y_POS_rot/Rs,z_POS_rot/Rs,'*'); grid on; xlabel('x'); hold on; ylabel('y'); zlabel('z'); axis equal;
%  %       plot3(x_obs_rot/Rs,y_obs_rot/Rs,z_obs_rot/Rs,'*'); %xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
%            plot3(r_LOS_steps(1,:)/Rs,r_LOS_steps(2,:)/Rs,r_LOS_steps(3,:)/Rs,'color','g'); hold on;
%       
% 
%            plot3(x_data/Rs,y_data/Rs,z_data/Rs,'.','color','k');
%             plot3(x_data(grid_indices(ii,isfinite(grid_indices(ii,:))))/Rs,y_data(grid_indices(ii,isfinite(grid_indices(ii,:))))/Rs,z_data(grid_indices(ii,isfinite(grid_indices(ii,:))))/Rs,'o','color','b');
%         
%             pause; close;
%             %
%         %     keyboard;
%         %
%         %     %%
                
        % Coordinates at which G-factors are evaluated.
        x_Gfactor = x_data(grid_indices(ii,isfinite(grid_indices(ii,:))));
        y_Gfactor = y_data(grid_indices(ii,isfinite(grid_indices(ii,:))));
        z_Gfactor = z_data(grid_indices(ii,isfinite(grid_indices(ii,:))));
        
        % Thomson scattering parameters along the LOS at the grid points.
        [G_T_tmp,G_P_tmp,G_R_tmp,G_tot_tmp] = G_Thomson(x_obs_rot,y_obs_rot,z_obs_rot,x_Gfactor,y_Gfactor,z_Gfactor,u);
        
        G_T_LOS(ii,1:length(G_T_tmp)) = G_T_tmp;
        G_P_LOS(ii,1:length(G_P_tmp)) = G_P_tmp;
        G_R_LOS(ii,1:length(G_R_tmp)) = G_R_tmp;
        G_tot_LOS(ii,1:length(G_tot_tmp)) = G_tot_tmp;
        
    end;
    
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

