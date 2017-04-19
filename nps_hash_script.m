% hash_script
% temp hack to coordinate GCS model cube with correct meshgrid

% NPS 10 Apr 2017

rmin = 0; % 12 used by AP
rmax = 15; % 12 used by AP

nElements = 128;
%  dGrid = 0.125; % Rs. rough desire
dGrid = (rmax - rmin)/(nElements-1);

[x_rs,y_rs,z_rs] = meshgrid(rmin:dGrid:rmax,rmin:dGrid:rmax,rmin:dGrid:rmax); 
[n,m,k] = size(x_rs);

x_data = x_rs * rs2km *1000; 
y_data = y_rs * rs2km *1000; 
z_data = z_rs * rs2km *1000; 
clear x_rs y_rs z_rs

load('basic_shape_plusnoise128.mat')
data = O_Struc.D1;

save('ProxyCubeDataTest','data','x_data','y_data','z_data')

% Tomography_Using_ART

% Vq = interp3(X,Y,Z,data_gcs,Xq,Yq,Zq);




