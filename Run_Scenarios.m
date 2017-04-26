%RUN_SCENARIOS. Run different tomographic reconstruction scenarios.

cc; tic;

% 1 AU
AU = 149598000e3; % m.

%% Main settings.

% Load the grid and true electron density.
load DataCube_CurrentSheet;
% Scenario save directory path and name.
scenario.name = './CurrentSheet_run';

% Coronagraph field-of-view and resolution.
%scenario.FOV = [2 17]; % Rs. Spherical & Gaussian.
scenario.FOV = [2 9]; % Rs.
scenario.resolution = 50; % arcsec.

% Heliocentric radius of the observations.
scenario.r_obs = 1*AU; % m.
% Heliocentric longitude of the observations.
scenario.theta_obs = linspace(-70,70,6); % deg.
scenario.lambda_obs = zeros(size(scenario.theta_obs)); % deg. Heliocentric latitude.

% The number of iterations.
scenario.no_of_ART_iterations = 3;

%%

% True data and reconstructions grid.
scenario.x = x_data; % m.
scenario.y = y_data; % m.
scenario.z = z_data; % m.
% The true electron density.
scenario.electron_density = data; % #/m^3.

% Limb darkening coefficient.
scenario.limb_darkening = 0.56;

% Thomson G-factor to be used.
scenario.G_factor = 'G_tot';

% Plot the output from the scenario.
scenario.plot = 1;

% Slices used in the plotting of the 3D volume.
scenario.slices_plot = slices_plot;

% Coordinates of the line plot through the cube.
scenario.line_plot = line_plot;

% Run the tomography script.
run_Tomography_Scenarios(scenario);

toc;
