%TOMOGRAPHY_USING_ART Try using algebraic reconstruction technique.

tic

cc;

%%

% Solar radius.
Rs = 695700e3; % m.
% 1 AU
AU = 149598000e3; % m.

FOV = [2 15]; % Rs.
resolution = 200; % arcsec.

r_obs = 1*AU; % m.
theta_obs =  linspace(-70,70,2); % deg. Heliocentric longitude. [-45 45] [-45 0 45] [-45 -15 15 45] [-45 -20 -10 10 20 45] [-60 -45 -20 -10 10 20 45 60]
NcamViews = length(theta_obs); % Number of camera view points. reduces single calc thorugh for loops
lambda_obs = zeros(size(theta_obs)); % deg. Heliocentric latitude.

% The number of iterations.
no_of_ART_iterations = 1;

% Limb darkening coefficient.
u = 0.56;
%
load CubeDataTest; % load ProxyCubeDataTest;

% Thomson G-factor to be used.
G_factor = 'G_tot';


%%
% Data dimension.
[i_dim,j_dim,k_dim] = size(x_data);

data_LOS_master = [];
G_tot_master = [];
G_tot_tmp_master = zeros(1,length(x_data(:)));

rms_difference = zeros(1,no_of_ART_iterations*NcamViews);
analysis_counter = 1;

% Initialize the reconstructed data cube. Use the same grid as for the
% "true" data.
Ne_inverted = zeros(size(data));



% %% CALC n record brensenham indices FOR ALL PIXELS and FOR ALL CAMERAS.
% 
% for ll = 1: NcamViews
%     for pp = 1:N_pix_in_image
%         % some more setup stuff needed
%         imageInd.OUT_Indices{1,ll}(:,pp) = NPS_eg_bresen (Ist_loc, Ien_loc,IN_extra);
%     end
% end
% 


%% multiple iterate through all camera view locations
for art_iterations = 1:no_of_ART_iterations
    fprintf('   ART iteration %01.0f/%01.0f...\n',art_iterations,no_of_ART_iterations);
    %% single Iterate Reconstruction over the viewpoints.
    for ii = 1:NcamViews
        fprintf('     Processing observation %01.0f/%01.0f...\n',ii,NcamViews);
        
        % Container for treated grid indices in the reconstruction cube.
        covered_grid_indices = zeros(1,length(x_data(:))); grid_container_counter = 1;
        
        % Generate synthetic coronagraph image from "true" data.
        [y_POS,z_POS,resolution_meters,data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,data,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
                
        % Generate coronagraph image from data volume being populated. This
        % step could be skipped for the first iteration but is included here
        % for simplicity and laziness.
        [y_POS_iteration,z_POS_iteration,resolution_meters_iteration,data_2D_LOS_iteration] = generate_2D_LOS_data(x_data,y_data,z_data,Ne_inverted,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
        
        % Compute the G factors and the LOS points within the grid.
        [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = map_LOS_2_G_data(x_data,y_data,z_data,y_POS,z_POS,r_obs,theta_obs(ii),lambda_obs(ii),u);
        
        % Loop over the LOS to assign the new values to the grid points.
        for ss = 1:length(y_POS),
            
            % Associated grid points in a given LOS.
            tmp_grid_indices = grid_indices(ss,:);
            % Finite indices. This gets a bit convoluted as we have indices for
            % the grid and indices for the matrix that holds the grid indices.
            finite_indices = isfinite(tmp_grid_indices);
            finite_grid_indices = tmp_grid_indices(finite_indices);
            
            % Use only grid points that were not populated by another LOS. Some
            % convolution with the indices here per the comment above...
            %[finite_grid_indices,finite_indices] = setdiff(finite_grid_indices,covered_grid_indices);
            [finite_grid_indices,finite_indices] = LIGHT_setdiff(finite_grid_indices,covered_grid_indices);
            
            % Continue only of the LOS "pierced" the data cube at more than one data point.
            if length(finite_grid_indices) > 1,
                
                % Add the new indices to the treated grid indices container.
                covered_grid_indices(grid_container_counter:grid_container_counter + length(finite_grid_indices) - 1) = finite_grid_indices;
                % Add to the counter.
                grid_container_counter = grid_container_counter + length(finite_grid_indices);
                
                % Segment lengths of the grid points along the LOS ray. YOU MAY WANT TO CHECK THIS MORE!
                dr_LOS_ray = cube_pierce_length(ss)/(length(finite_grid_indices) - 1);
                % Distribution of the electrons per ART rule. Weights based on Thomson parameters and the segment lengths.
                data_spread_along_LOS = (data_2D_LOS(ss) - data_2D_LOS_iteration(ss))*(1/length(finite_grid_indices)).*(1./G_tot_LOS(ss,finite_indices))*(1/dr_LOS_ray);
                Ne_inverted(finite_grid_indices) = Ne_inverted(finite_grid_indices) + data_spread_along_LOS;
                                
            end;
            
            
        end;
        
        % Smoothing of the data.
        Ne_inverted = ordfilt3D(Ne_inverted,14);
        
        % TEST TEST
        % Ne_inverted_gridded = reshape(Ne_inverted,31,31,31);
        Ne_inverted_ART_gridded = reshape(Ne_inverted,i_dim,j_dim,k_dim);
        
        figure; slice(x_data/Rs,y_data/Rs,z_data/Rs,Ne_inverted_ART_gridded,[0],[5],[-3.5]); colorbar; title('Reconstructed electron density [#/m^3] in the solar corona'); xlabel('x [Rs]'); ylabel('y [Rs]'); zlabel('z [Rs]');
        caxis([0 max(data(:))]); xlim([min(x_data(:)) max(x_data(:))]/Rs); ylim([min(y_data(:)) max(y_data(:))]/Rs); zlim([min(z_data(:)) max(z_data(:))]/Rs);
        colormap(gray); shading('interp');
        
        [y_gridded, z_gridded] = meshgrid(min(y_POS):resolution_meters:max(y_POS),min(z_POS):resolution_meters:max(z_POS));
         
        data_2D_LOS_grid = griddata(y_POS,z_POS,data_2D_LOS,y_gridded,z_gridded,'linear');
        r_gridded = sqrt(y_gridded.^2 + z_gridded.^2);
        kk = find(r_gridded < FOV(1)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
        for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
        
        figure; pcolor(y_gridded/Rs,z_gridded/Rs,data_2D_LOS_grid); colorbar; caxis([min(data_2D_LOS_grid(:)) max(data_2D_LOS_grid(:))]);
        title(sprintf('Iteration at %01.0f helio longitude, %01.0f latitude',theta_obs(ii),lambda_obs(ii)));
        xlabel('y [Rs]'); ylabel('z [Rs]');
        colormap(gray); shading('interp');
        
        %   eval(sprintf('print -dpng ./plots/2D_LOS_data_%01.0f.png',ii)); close;
        
        % RMS difference.
        rms_difference(analysis_counter) = sum( sqrt( ( data(:) - Ne_inverted_ART_gridded(:) ).^2 ) )/length(data(:));
        analysis_counter = analysis_counter + 1;
        
        disp(sprintf('     %01.0f/%01.0f spacecraft. RMS difference in the electron density %01.0f #/m^3.',ii,NcamViews,rms_difference(analysis_counter-1)));
        
        
    end; % Different viewpoints.
    
end; % ART iterations.

figure; slice(x_data/Rs,y_data/Rs,z_data/Rs,data,[0],[5],[-3.5]); colorbar; title('"True" electron density [#/m^3] in the solar corona'); xlabel('x [Rs]'); ylabel('y [Rs]'); zlabel('z [Rs]');
caxis([0 max(data(:))]); xlim([min(x_data(:)) max(x_data(:))]/Rs); ylim([min(y_data(:)) max(y_data(:))]/Rs); zlim([min(z_data(:)) max(z_data(:))]/Rs);
colormap(gray); shading('interp');

figure; plot(squeeze(z_data(26,26,:)/Rs),squeeze(Ne_inverted(26,26,:)),'k'); hold on; plot(squeeze(z_data(26,26,:)/Rs),squeeze(data(26,26,:)),'b');
xlabel('z [Rs]'); ylabel('Electron density [#/m^3]'); title('True (blue) and reconstructed (black) electron density'); grid on;

eval(sprintf('print -dpng ./plots/Line_plot_true_reconstructed_%01.0f_spacecraft.png',NcamViews)); close;

figure; plot(rms_difference); grid on;
title(sprintf('Convergence for %1.0f ART iterations and %01.0f spacecraft',no_of_ART_iterations,NcamViews)); xlabel('ART step'); ylabel('RMS [#/m^3]');

%
% figure; slice(x_data/Rs,y_data/Rs,z_data/Rs,Ne_inverted_svd_gridded,[4],[5],[5]); colorbar; title('Reconstructed electron density [#/m^3] in the solar corona'); xlabel('x [Rs]'); ylabel('y [Rs]'); zlabel('z [Rs]');
% caxis([0 max(data(:))]); xlim([min(x_data(:)) max(x_data(:))]/Rs); ylim([min(y_data(:)) max(y_data(:))]/Rs); zlim([min(z_data(:)) max(z_data(:))]/Rs);
% colormap(gray); shading('interp');
%
% for ii = 1:length(theta_obs)
%
%     % Generate coronagraph images from the reconstructed data.
%     [y_POS,z_POS,resolution_meters,data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,Ne_inverted_svd_gridded,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,dr,u,G_factor);
%
%     data_2D_LOS_grid = griddata(y_POS,z_POS,data_2D_LOS,y_gridded,z_gridded,'linear');
%     r_gridded = sqrt(y_gridded.^2 + z_gridded.^2);
%     kk = find(r_gridded < FOV(1)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
%     for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
%
%     figure; pcolor(y_gridded/Rs,z_gridded/Rs,data_2D_LOS_grid); colorbar; caxis([0 5e-12]);
%     title(sprintf('Reconstruction at %01.0f helio longitude, %01.0f latitude',theta_obs(ii),lambda_obs(ii)));
%     xlabel('y [Rs]'); ylabel('z [Rs]');
%     colormap(gray); shading('interp');
%
%     eval(sprintf('print -dpng ./plots/2D_reconstructed_data_%01.0f.png',ii)); close;
%
%
% end;

toc
