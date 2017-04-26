%CONSTRUCT_3D_DATA_CUBE
%

%   Antti Pulkkinen, March 2017.

cc;

% The structure type;
structureType = 'Spherical&Gaussian';         % CurrentSheet   Gaussian

% Solar radius.
Rs = 695700e3; % m.
% 1 AU
AU = 149598000e3; % m.

% Grid definition.
dGrid = 0.5; % Rs.
nElements = 40;

% The grid;
%[x,y,z] = meshgrid(-dGrid*nElements/2:dGrid:dGrid*nElements/2,-dGrid*nElements/2:dGrid:dGrid*nElements/2,-dGrid*nElements/2:dGrid:dGrid*nElements/2);
% Initialize data.
%Ne = zeros(size(x));

switch structureType
    
    case 'Gaussian'
               
        [x,y,z] = meshgrid(0:0.5:12,0:0.5:12,0:0.5:12); [n,m,k] = size(x);
        Ne = zeros(size(x));       
        
        %Ne = 11^10*exp(-(1e-1*(x - 5).^2 + 5e-1*(y - 10).^2 + 1e-1*(z - 5).^2)) + 10^10*exp(-(1e-1*(x + 5).^2 + 1e-1*(y - 5).^2 + 1e-1*(z - 5).^2)) ;
        Ne = 11^10*exp(-(2e-1*(x + 2).^2 + 25e-1*(y - 5).^2 + 1e-1*(z - 5).^2)) + 10^10*exp(-(1e-1*(x + 5).^2 + 1e-1*(y - 4).^2 + 1e-1*(z - 5).^2)) ;
        %Ne = 11^10*exp(-(2e-1*(x + 2).^2 + 25e-1*(y - 5).^2 + 1e-1*(z - 5).^2)) - 10^10*exp(-(1e-1*(x + 5).^2 + 1e-1*(y - 4).^2 + 1e-1*(z - 5).^2)) ;
        
    case 'Spherical'
        
        
        [x,y,z] = meshgrid(0:0.5:12,0:0.5:12,0:0.5:12); [n,m,k] = size(x);
        Ne = zeros(size(x));        
        
        Ne = 10^9*((x - 4).^2 + (y - 4).^2 + (z - 4).^2);
        sphere_center = [6 6 6]; sphere_radius = 5;
        radius_from_sphere_center = sqrt( (x - sphere_center(1)).^2 + (y - sphere_center(2)).^2 + (z - sphere_center(3)).^2 );
        kk = find(radius_from_sphere_center > sphere_radius);
        Ne(kk) = 0; Ne = reshape(Ne,n,m,k);
        
    case 'Spherical&Gaussian'
        
        [x,y,z] = meshgrid(0:0.15:12,0:0.15:12,0:0.15:12); [n,m,k] = size(x);
        Ne = zeros(size(x));
        
        Ne = 10^9*((x - 4).^2 + (y - 4).^2 + (z - 4).^2);
        sphere_center = [6 6 6]; sphere_radius = 5;
        radius_from_sphere_center = sqrt( (x - sphere_center(1)).^2 + (y - sphere_center(2)).^2 + (z - sphere_center(3)).^2 );
        kk = find(radius_from_sphere_center > sphere_radius);
        Ne(kk) = 0; Ne = reshape(Ne,n,m,k);
        % Add Gaussian small-scale structure.
        Ne_gauss = 11^10*exp(-(20e-1*(x - 5).^2 + 20e-1*(y - 5).^2 + 20e-1*(z - 7).^2));
        Ne = Ne + Ne_gauss;
        
        % Slices used in the plotting.
        slices_plot = {5, 5, 7};
        % Coordinates of the line plot through the cube.
        line_plot = {'34,34,:'};
        
    case 'CurrentSheet'
        
        % Grid resolution.
        dr_grid = 0.1; % Rs.
        
        [x,y,z] = meshgrid(2:dr_grid:8,-3:dr_grid:3,-3:dr_grid:3); [n,m,k] = size(x);
        Ne = zeros(size(x));
        
        % Sheet.
        kk_x = find(x <= 7 & x >= 3.5); kk_y = find(y <= 1 & y >= -1); kk_z = find(z <= 0.2 & z >= -0.2);
        kk_tmp = intersect(kk_x,kk_y); kk_all = intersect(kk_tmp,kk_z);
        Ne(kk_all) = 10^10;
        
        % Add speherical blob(s).
        radius_blob_1 = 0.6;
        r_blob_1 = sqrt( (x - 7).^2 + (y - 0).^2 + (z - 0).^2 );
        kk_blob_1 = find (r_blob_1 <= radius_blob_1);
        Ne(kk_blob_1) = 1.3*(10^10);
        
        radius_blob_2 = 0.4;
        r_blob_2 = sqrt( (x - 5).^2 + (y - 0).^2 + (z - 0).^2 );
        kk_blob_2 = find (r_blob_2 <= radius_blob_2);
        Ne(kk_blob_2) = 1.3*(10^10);
        
        % Rotate the entire structure. Without rotations small scale
        % structure in ecliptic plane will be obstructed by occulter from
        % many viewpoints.
        
        % Rotation matrices for rotation about z- and y-axis (theta and lambda, respectively).
        theta = 90; % longitude rotation.
        lambda = -45; % latitude rotation.
        
        theta_rad = theta*pi/180; lambda_rad = lambda*pi/180;
        R_z = [cos(theta_rad) -sin(theta_rad) 0 ; sin(theta_rad) cos(theta_rad) 0 ; 0 0 1];
        R_y = [cos(-lambda_rad) 0 sin(-lambda_rad) ; 0 1 0 ; -sin(-lambda_rad) 0 cos(-lambda_rad)];
        
        % Rotate the plane of sky to correspond to the location of the observer.
        for ii = 1:length(x(:)),
            
            r_rot_zy = R_z*R_y*[x(ii) ; y(ii) ; z(ii)];
            
            x(ii) = r_rot_zy(1);
            y(ii) = r_rot_zy(2);
            z(ii) = r_rot_zy(3);
            
        end;
        
        % Regenerate the grid for the structure to be proper meshgrid.
        [x_new,y_new,z_new] = meshgrid(min(x(:)):dr_grid:max(x(:)),min(y(:)):dr_grid:max(y(:)),min(z(:)):dr_grid:max(z(:)));
        % Interpolate.
        Ne_new = griddata(x,y,z,Ne,x_new,y_new,z_new,'nearest');
        % Replace the original variable.
        Ne = Ne_new;
        x = x_new; y = y_new; z = z_new;
        
end;

figure; slice(x,y,z,Ne,slices_plot{:}); colorbar; colormap('gray'); shading('interp'); xlabel('x [Rs]'); ylabel('y [Rs]'); zlabel('z [Rs]'); xlim([-15 15]); ylim([-15 15]); zlim([-15 15]);

figure; plot(eval(sprintf('squeeze(z(%s))/Rs',line_plot{:})),eval(sprintf('squeeze(Ne(%s))',line_plot{:})),'k');
xlabel('distance [Rs]'); ylabel('Electron density [#/m^3]'); title('Electron density through the volume'); grid on;

x_data = x*Rs; y_data = y*Rs; z_data = z*Rs; data = Ne; save CubeDataTest x_data y_data z_data data;


