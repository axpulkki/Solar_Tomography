%
%   UNITS of operation in Rs
%
% test script. simply run routine in empty workspace

% test_NPS_eg_bresen

% ray path preamble - ie origins of camera pixel
orig_r = 0.1 * au2rs ; % * rs2km *1000;
y_r = 4; z_r =1;

%% create camera position ie ray origin/end
orig = [orig_r, y_r, z_r];

Ist_loc = orig;
enloc = -1 * orig;       % reflection about the Sun 
Ien_loc = enloc - ([z_r, 0.5*y_r, z_r]);

%% Define cube size/location
load CubeDataGaussian
IN_extra.cube = data;
IN_extra.xgrid = x_data;        % x = x_data/695700e3;
IN_extra.ygrid = y_data;
IN_extra.zgrid = z_data;

%% fid bresenham indices through cube
OUT_Indices = NPS_eg_bresen (Ist_loc, Ien_loc, IN_extra);


