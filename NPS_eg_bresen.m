function OUT_Indices = NPS_eg_bresen (Ist_loc, Ien_loc,IN_extra)
% 
%   Ist_loc     =   Start position of single ray. ie camera pixel 
%                   location in phyisical space.
% 
%   Ien_loc     =   End position of single ray. ie Arbitrary position
%                   infinitely far away along correct ray direction
% 
%   IN_extra    = [Structure]. values to incl. any other required inputs
%                   xgrid, ygrid, zgrid of data cube
% 
%   OUT_Indices = index array of ea voxel for single ray path though data cube
%
%   UNITS of operation in Rs

% try with test_NPS_eg_bresen
% NPS. 17 Apr 2017 


st_loc = Ist_loc ; % /1000) * km2rs;  
en_loc = Ien_loc ; % /1000) * km2rs;

n_f = 10;
orig = st_loc;              % ray's origin, phys. space
dir = (en_loc - orig)/n_f;      % ray's direction, phys. space

stex = IN_extra;
xg = (stex.xgrid/1000)* km2rs; yg = (stex.ygrid/1000)* km2rs; zg = (stex.zgrid/1000)* km2rs;
cube = stex.cube;



%%%%%%%% OBJECTIVE IS TO FIND WHERE RAY %%%%%%%%%%%%%%%
%%%%%%%%    INTERSECTS WITH DATA CUBE, Physical space   %%%%%%%%%%%%%%%

%% Find all apex of data cube in physical space.
x(1)= min(min(min(xg))); x(2) = max(max(max(xg)));
y(1)= min(min(min(yg))); y(2) = max(max(max(yg)));
z(1)= min(min(min(zg))); z(2) = max(max(max(zg)));

% apex_all = [x(1), x(1), x(2), x(2), x(1), x(2), x(2), x(1);...
%             y(1), y(2), y(2), y(1), y(1), y(2), y(2), y(1);...
%             z(1), z(1), z(1), z(1), z(2), z(2), z(2), z(2)];
% vert_all = apex_all';

% reduce mgrid
[xgd,ygd,zgd] = meshgrid(x(1):(x(2)-x(1)):x(2),y(1):(y(2)-y(1)):y(2),z(1):(z(2)-z(1)):z(2));
% [xgd,ygd,zgd] = meshgrid(x(1):(x(2)-x(1))/2:x(2),y(1):(y(2)-y(1))/2:y(2),z(1):(z(2)-z(1))/2:z(2));
DT = delaunayTriangulation(xgd(:), ygd(:), zgd(:));
% [xgd,ygd,zgd] = meshgrid(x,y,z);        [n,m,k] = size(xgd);
[faces, vertices] = freeBoundary(DT);
vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);

%% find ray intersect pts with ALL planes 
% (More than 6 planes of cube, as code uses traingular mesh)
[Intersect,Dist,bcd1,bcd2, Xcord] = TriangleRayIntersection(orig, dir, ...
                                            vert1, vert2, vert3);
f2 = figure;
trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),Intersect*1.0,'FaceAlpha', 0.7)
hold on

lend = n_f * dir;
line('XData',orig(1)+[0 lend(1)],'YData',orig(2)+[0 lend(2)],'ZData',...
  orig(3)+[0 lend(3)],'Color','r','LineWidth',3)
axis equal
% set(gca, 'CameraPosition', [106.2478  -35.9079  136.4875])

%% Find entry and exit to cube, in physical space
% ie find shortest and longest distance of interesects
inter_pos = Xcord(Intersect);
inter_dis = Dist(Intersect);

[entry_D,Ien_D] = min(inter_dis);
[exit_D ,Iex_D] = max(inter_dis);

entry_P = inter_pos(Ien_D,:);           % units still in Rs
exit_P  = inter_pos(Iex_D,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Find entry/ exit of cube as index cube value.






%% Find index array for ray along Data cube
[xpt,ypt,zpt] = bresenham_line3d (st,en,0);



OUT_Indices = 1;

return


