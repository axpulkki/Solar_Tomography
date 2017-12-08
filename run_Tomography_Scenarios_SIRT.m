function run_Tomography_Scenarios(scenario)
%run_Tomography_Scenarios_SIRT(scenario) ...
%

% Antti Pulkkinen, July 2017.

% Solar radius.
Rs = 695700e3; % m.

%% Scenario settings.

x_data = scenario.x; % m.
y_data = scenario.y; % m.
z_data = scenario.z; % m.
data = scenario.electron_density; % #/m^3.
if scenario.background_subtraction,
    data_background = scenario.electron_density_background; % #/m^3.
end;

FOV = scenario.FOV; % Rs.
resolution = scenario.resolution; % arcsec.

r_obs = scenario.r_obs; % m.
theta_obs =  scenario.theta_obs; % deg. Heliocentric longitude.
lambda_obs = scenario.lambda_obs; % deg. Heliocentric latitude.
NcamViews = length(theta_obs); % Number of camera view points.

% The number of iterations.
no_of_ART_iterations = scenario.no_of_ART_iterations;

% Limb darkening coefficient.
u = scenario.limb_darkening;

% Thomson G-factor to be used.
G_factor = scenario.G_factor;

%%

% % Data dimension.
% [i_dim,j_dim,k_dim] = size(x_data);

data_LOS_master = [];
G_tot_master = [];
G_tot_tmp_master = zeros(1,length(x_data(:)));

rms_difference = zeros(1,no_of_ART_iterations);
analysis_counter = 1;

% Initialize the reconstructed data cube. Use the same grid as for the true data.
Ne_inverted = zeros(size(data));

% Generate the scenario folder.
mkdir(scenario.name);

% Coordinates of the line and the distance from the first point in the line.
x_line = eval(sprintf('squeeze(x_data(%s))',scenario.line_plot{:}));
y_line = eval(sprintf('squeeze(y_data(%s))',scenario.line_plot{:}));
z_line = eval(sprintf('squeeze(z_data(%s))',scenario.line_plot{:}));
%line_plot_distance = sqrt( (x_line(1) - x_line).^2 + (y_line(1) - y_line).^2 + (z_line(1) - z_line).^2 );

% Generate the coronagraph images from the true data.
for ii = 1:NcamViews,
    
    fprintf('   RUN_TOMOGRAPHY_SCENARIOS: Generating coronagraph image %01.0f/%01.0f...\n',ii,NcamViews);
    
    % Generate synthetic coronagraph images from "true" data.
    [camView(ii).y_POS,camView(ii).z_POS,camView(ii).resolution_meters,camView(ii).data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,data,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
    
    % Images of the background.
    if scenario.background_subtraction,
        
        fprintf('   RUN_TOMOGRAPHY_SCENARIOS: Generating coronagraph image of the background %01.0f/%01.0f...\n',ii,NcamViews);
        [camView(ii).y_POS,camView(ii).z_POS,camView(ii).resolution_meters,camView(ii).data_2D_LOS_background] = generate_2D_LOS_data(x_data,y_data,z_data,data_background,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
        % Subtract the background from the image.
        camView(ii).data_2D_LOS = camView(ii).data_2D_LOS - camView(ii).data_2D_LOS_background;
        % Force the negative pixel values to zero.
        kk_LOS_negative = find( camView(ii).data_2D_LOS < 0 ); camView(ii).data_2D_LOS(kk_LOS_negative) = 0;
        
    end;
    
    if scenario.plot
        
        % Plot and save the coronagraph images.
        [y_gridded, z_gridded] = meshgrid(min(camView(ii).y_POS):camView(ii).resolution_meters:max(camView(ii).y_POS),min(camView(ii).z_POS):camView(ii).resolution_meters:max(camView(ii).z_POS));
        
        data_2D_LOS_grid = griddata(camView(ii).y_POS,camView(ii).z_POS,camView(ii).data_2D_LOS,y_gridded,z_gridded,'nearest');
        
        if scenario.background_subtraction,
            data_2D_LOS_grid_background = griddata(camView(ii).y_POS,camView(ii).z_POS,camView(ii).data_2D_LOS_background,y_gridded,z_gridded,'nearest');
        else
            data_2D_LOS_grid_background = 1;
        end;
        
        r_gridded = sqrt(y_gridded.^2 + z_gridded.^2);
        
        % Carve out the inner FOV.
        kk = find(r_gridded < FOV(1)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
        for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
        
        % Carve out the outer FOV.
        kk = find(r_gridded > FOV(2)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
        for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
        
        coronagraph_plot_quantity = data_2D_LOS_grid.*r_gridded.^2;
        
        figure; pcolor(y_gridded/Rs,z_gridded/Rs,coronagraph_plot_quantity); colorbar; caxis([min(coronagraph_plot_quantity(:)) max(coronagraph_plot_quantity(:))]);
        title(sprintf('Camera view at %01.0f deg. helio longitude, %01.0f deg. latitude.',theta_obs(ii),lambda_obs(ii)));
        xlabel('y [Rs]'); ylabel('z [Rs]');
        colormap(gray); shading('interp');
        
        eval(sprintf('print -dpng %s/coronagraph_image_%01.0f.png',scenario.name,ii)); close;
        
    end;
    
end;

% If background was subtracted from the analysis, compare below
% reconstructed electron densitites to true electron densities with
% background subtracted.
if scenario.background_subtraction,
    
    disp('    RUN_TOMOGRAPHY_SCENARIOS: subtracting background from the true electron densities...');
    data = data - data_background;
    % Set negative values to zero.
    data(find(data < 0)) = 0;
    
end;

%% multiple iterate through all camera view locations
for art_iterations = 1:no_of_ART_iterations
    fprintf('   RUN_TOMOGRAPHY_SCENARIOS: ART iteration %01.0f/%01.0f...\n',art_iterations,no_of_ART_iterations);
    
    % Initiate updates data cube.
    Ne_collect_updates = zeros(size(data));
    % Initiate the cube indicating the number of updates in each cell.
    Ne_number_of_updates = zeros(size(data));
    
    %% single Iterate Reconstruction over the viewpoints.
    for ii = 1:NcamViews
        
        fprintf('    RUN_TOMOGRAPHY_SCENARIOS: Processing observation %01.0f/%01.0f...\n',ii,NcamViews);
        
        % Container for treated grid indices in the reconstruction cube.
        covered_grid_indices = zeros(1,length(x_data(:))); grid_container_counter = 1;
        
        % Generate coronagraph image from data volume being populated. This
        % step could be skipped for the first iteration but is included here
        % for simplicity and laziness.
        [y_POS_iteration,z_POS_iteration,resolution_meters_iteration,data_2D_LOS_iteration] = generate_2D_LOS_data(x_data,y_data,z_data,Ne_inverted,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
                
        % Compute the G factors and the grid points intersected by the LOS ray.
        [G_T_LOS,G_P_LOS,G_R_LOS,G_tot_LOS,grid_indices,cube_pierce_length] = map_LOS_2_G_data(x_data,y_data,z_data,camView(ii).y_POS,camView(ii).z_POS,r_obs,theta_obs(ii),lambda_obs(ii),u);
        
        % Loop over the LOS to assign the new values to the grid points.
        for ss = 1:length(camView(ii).y_POS),
            
            % Associated grid points in a given LOS.
            tmp_grid_indices = grid_indices(ss,:);
            % Finite indices. This gets a bit convoluted as we have indices for
            % the grid and indices for the matrix that holds the grid indices.
            finite_indices = isfinite(tmp_grid_indices);
            finite_grid_indices = tmp_grid_indices(finite_indices);
            
%             % Use only grid points that were not populated by another LOS. Some
%             % convolution with the indices here per the comment above...
%             %[finite_grid_indices,finite_indices] = setdiff(finite_grid_indices,covered_grid_indices);
%             [finite_grid_indices,finite_indices] = LIGHT_setdiff(finite_grid_indices,covered_grid_indices);
            
            % Continue only of the LOS "pierced" the data cube at more than one data point.
            if length(finite_grid_indices) > 1,
                
                
%                 % Add the new indices to the treated grid indices container.
%                 covered_grid_indices(grid_container_counter:grid_container_counter + length(finite_grid_indices) - 1) = finite_grid_indices;
%                 % Add to the counter.
%                 grid_container_counter = grid_container_counter + length(finite_grid_indices);
                                
                % Segment lengths of the grid points along the LOS ray. WE NEED TO CHECK THAT THIS IS CORRECT!
                dr_LOS_ray = cube_pierce_length(ss)/(length(finite_grid_indices) - 1);
                % Distribution of the electrons per ART rule. Weights based on Thomson parameters and the segment lengths.
                data_spread_along_LOS = (camView(ii).data_2D_LOS(ss) - data_2D_LOS_iteration(ss))*(1/length(finite_grid_indices)).*(1./G_tot_LOS(ss,finite_indices))*(1/dr_LOS_ray);
                
                % Hamming window.
                N = length(data_spread_along_LOS);
                hamming_window = 0.54 - 0.46*cos( 2*pi*(0:(N-1))/(N-1) );               
                
                % Collect updates. Apply Hamming window. Note that multiple LOS can pierce a given cell.
                Ne_collect_updates(finite_grid_indices) = Ne_collect_updates(finite_grid_indices) + data_spread_along_LOS.*hamming_window;
                % Update the number of updates of each cell.
                Ne_number_of_updates(finite_grid_indices) = Ne_number_of_updates(finite_grid_indices) + 1;
                
            end;
            
        end;
        
        % Find the cells that were updated.
        indices_cells_updated = find(Ne_number_of_updates > 0);
        
        % SIRT update.
        Ne_inverted(indices_cells_updated) = Ne_inverted(indices_cells_updated) + Ne_collect_updates(indices_cells_updated)./Ne_number_of_updates(indices_cells_updated);
           
        % Force positivity of the result. Set negative values to zero.
        kk_negative = find(Ne_inverted < 0); Ne_inverted(kk_negative) = 0;

    end; % Different viewpoints.
    
    % Smoothing of the data. WE NEED TO TEST ALSO OTHER SMOOTHING ALGORITHMS!
    Ne_inverted = ordfilt3D(Ne_inverted,14);
    
    
    % RMS difference.
    rms_difference(analysis_counter) = sum( sqrt( ( data(:) - Ne_inverted(:) ).^2 ) )/length(data(:));
    analysis_counter = analysis_counter + 1;
    
    fprintf('    RUN_TOMOGRAPHY_SCENARIOS: After iteration %01.0f/%01.0f RMS difference in the electron density %01.0f #/m^3.\n',art_iterations,no_of_ART_iterations,rms_difference(analysis_counter-1));

    
end; % ART iterations.

% Save the run data.
eval(sprintf('save %s/run_data.mat',scenario.name));

if scenario.plot
    
    [n,m,s] = size(data);
    
    eval(sprintf('mkdir %s/slices_frames',scenario.name'));
    
    for ssSlices = 1:s,
        
        figure; 
        subplot(1,2,1), pcolor(squeeze(x_data(:,:,ssSlices))/Rs,squeeze(y_data(:,:,ssSlices))/Rs,(squeeze(data(:,:,ssSlices)))); 
        colorbar(gca,'SouthOutside'); title(sprintf('True z: %01.1f Rs',z_data(1,1,ssSlices)/Rs)); axis equal;
        caxis([min(min(squeeze(data(:,:,ssSlices)))) max(max(squeeze(data(:,:,ssSlices))))]);
        xlabel('x [Rs]'); ylabel('y [Rs]');
        
        subplot(1,2,2), pcolor(squeeze(x_data(:,:,ssSlices))/Rs,squeeze(y_data(:,:,ssSlices))/Rs,(squeeze(Ne_inverted(:,:,ssSlices)))); 
        colorbar(gca,'SouthOutside'); title(sprintf('Reconstructed z: %01.1f Rs',z_data(1,1,ssSlices)/Rs)); axis equal;
        caxis([min(min(squeeze(data(:,:,ssSlices)))) max(max(squeeze(data(:,:,ssSlices))))]);
        xlabel('x [Rs]'); ylabel('y [Rs]');
        
        eval(sprintf('print -dpng ./%s/slices_frames/frame_%01.0f.png',scenario.name,ssSlices)); close;
        
    end;
    
    figure; slice(x_data/Rs,y_data/Rs,z_data/Rs,log10(data),scenario.slices_plot{:}); colorbar; title('True electron density [#/m^3].'); xlabel('x [Rs]'); ylabel('y [Rs]'); zlabel('z [Rs]');
    %caxis([0 max(data(:))]); 
    xlim([min(x_data(:)) max(x_data(:))]/Rs); ylim([min(y_data(:)) max(y_data(:))]/Rs); zlim([min(z_data(:)) max(z_data(:))]/Rs);
    colormap(gray); shading('interp'); hold on;
    plot3(x_line/Rs,y_line/Rs,z_line/Rs,'linewidth',3);
    
    eval(sprintf('print -dpng %s/true_electron_density.png',scenario.name));
    
    figure; plot(rms_difference); grid on;
    title(sprintf('Convergence for %1.0f ART iterations and %01.0f spacecraft',no_of_ART_iterations,NcamViews)); xlabel('ART step'); ylabel('RMS [#/m^3]');
    
    eval(sprintf('print -dpng ./%s/ART_convergence.png',scenario.name));
    
    % Plot the coronagraph images obtained from the reconstructed electron density.
    for ii = 1:NcamViews,
        
        fprintf('   RUN_TOMOGRAPHY_SCENARIOS: Generating coronagraph image from the final reconstucted data %01.0f/%01.0f...\n',ii,NcamViews);
        
        % Generate synthetic coronagraph images from the reconstructed data.
        [camView(ii).y_POS,camView(ii).z_POS,camView(ii).resolution_meters,camView(ii).data_2D_LOS] = generate_2D_LOS_data(x_data,y_data,z_data,Ne_inverted,r_obs,theta_obs(ii),lambda_obs(ii),FOV,resolution,u,G_factor);
        
        % Plot and save the coronagraph images.
        [y_gridded, z_gridded] = meshgrid(min(camView(ii).y_POS):camView(ii).resolution_meters:max(camView(ii).y_POS),min(camView(ii).z_POS):camView(ii).resolution_meters:max(camView(ii).z_POS));
        
        data_2D_LOS_grid = griddata(camView(ii).y_POS,camView(ii).z_POS,camView(ii).data_2D_LOS,y_gridded,z_gridded,'nearest');
        r_gridded = sqrt(y_gridded.^2 + z_gridded.^2);
        
        % Carve out the inner FOV.
        kk = find(r_gridded < FOV(1)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
        for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
        
        % Carve out the outer FOV.
        kk = find(r_gridded > FOV(2)*Rs); [kk_i,kk_j] = ind2sub(size(r_gridded),kk);
        for jj = 1:length(kk), data_2D_LOS_grid(kk_i(jj),kk_j(jj)) = NaN; end;
        
        coronagraph_plot_quantity = data_2D_LOS_grid.*r_gridded.^2;
        
        figure; pcolor(y_gridded/Rs,z_gridded/Rs,coronagraph_plot_quantity); colorbar; caxis([min(coronagraph_plot_quantity(:)) max(coronagraph_plot_quantity(:))]);
        title(sprintf('Reconstr. camera view at %01.0f deg. helio longitude, %01.0f deg. latitude.',theta_obs(ii),lambda_obs(ii)));
        xlabel('y [Rs]'); ylabel('z [Rs]');
        colormap(gray); shading('interp');
        
        eval(sprintf('print -dpng %s/coronagraph_image_from_reconstruction_%01.0f.png',scenario.name,ii)); close;
        
    end;
    
end;


