%RUN_SCENARIOS. Run different tomographic reconstruction scenarios.

% Solar radius.
Rs = 695700e3; % m.
% 1 AU
AU = 149598000e3; % m.

% Load the grid and true electron density.
load CubeDataTest;

% Scenario save directory path and name.
scenario.name = './testScenario';

% True data and reconstructions grid.
scenario.x = x_data; % m.
scenario.y = y_data; % m.
scenario.z = z_data; % m.
% The true electron density.
scenario.electron_density = data; % #/m^3.

% Coronagraph field-of-view and resolution.
scenario.FOV = [2 17]; % Rs.
scenario.resolution = 150; % arcsec.

% Heliocentric radius of the observations.
scenario.r_obs = 1*AU; % m.
% Heliocentric longitude of the observations.
scenario.theta_obs = linspace(-70,70,6); % deg.
scenario.lambda_obs = zeros(size(scenario.theta_obs)); % deg. Heliocentric latitude.

% The number of iterations.
scenario.no_of_ART_iterations = 3;

% Limb darkening coefficient.
scenario.limb_darkening = 0.56;

% Thomson G-factor to be used.
scenario.G_factor = 'G_tot';

% Plot the output from the scenario.
scenario.plot = 1;

% Slices used in the plotting of the 3D volume.
scenario.slices_plot = {5, 5, 7};

% Coordinates of the line plot through the cube.
scenario.line_plot = {'34,34,:'};

% Run the tomography script.
run_Tomography_Scenarios(scenario);

