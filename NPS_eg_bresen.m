function OUT_Indices = NPS_eg_bresen (Ist_loc, Ien_loc,IN_extra)
% 
%   Ist_loc     =   Start position of single ray. ie camera pixel 
%                   location in phyisical space. [x,y,z] in units of [m]
% 
%   Ien_loc     =   End position of single ray. eg Arbitrary position
%                   infinitely far away along correct ray direction
%                   [x,y,z] in units of [m]
% 
%   IN_extra    = [Structure]. values to incl. any other required inputs
%                  =  xgrid, ygrid, zgrid of data cube (coming from meshgrid)
%                  = debug = [1] or [0] for figure plotting. default = 0
% 
%   OUT_Indices = index array of ea voxel for single ray path though data cube
%
%   UNITS of IN/OUT-puts [m]  -- code reconverts to Rs for ease of debugging 

% try with test_NPS_eg_bresen
% NPS. 17 Apr 2017 


st_loc = Ist_loc *(km2rs/1000); % /1000) * km2rs;  
en_loc = Ien_loc *(km2rs/1000); % /1000) * km2rs;

n_f = 10;
orig = st_loc;              % ray's origin, phys. space
dir = (en_loc - orig)/n_f;      % ray's direction, phys. space

stex = IN_extra;
xg = stex.xgrid; yg = stex.ygrid; zg = stex.zgrid;

xgdpts = xg(1,:,1) *(km2rs/1000) ;
ygdpts = yg(:,1,1)' * (km2rs/1000);
zgdpts = squeeze(zg(1,1,:))' * (km2rs/1000);
% [xg,yg,zg] = meshgrid(0:0.5:12,0:1:24,0:2:48); 
% xgdpts = xg(1,:,1); ygdpts = yg(:,1,1)'; zgdpts = squeeze(zg(1,1,:))';
% temp = xg(1,:,1);temp = yg(:,1,1)';temp = squeeze(zg(1,1,:))';

Fplot = stex.debug;
if (Fplot ~=0 & Fplot~=1); Fplot = 0; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OBJECTIVE IS TO FIND WHERE RAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    INTERSECTS WITH DATA CUBE, Physical space  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all apex of data cube in physical space.
% x(1)= min(min(min(xg))); x(2) = max(max(max(xg)));   y(1)= min(min(min(yg))); y(2) = max(max(max(yg)));   z(1)= min(min(min(zg))); z(2) = max(max(max(zg)));
x(1)= xgdpts(1); x(2) = xgdpts(end);
y(1)= ygdpts(1); y(2) = ygdpts(end);
z(1)= zgdpts(1); z(2) = zgdpts(end);    

%% Find vertices and faces of cube
% reduce mgrid
[xgd,ygd,zgd] = meshgrid(x(1):(x(2)-x(1)):x(2),y(1):(y(2)-y(1)):y(2),z(1):(z(2)-z(1)):z(2));
% [xgd,ygd,zgd] = meshgrid(x(1):(x(2)-x(1))/2:x(2),y(1):(y(2)-y(1))/2:y(2),z(1):(z(2)-z(1))/2:z(2));
DT = delaunayTriangulation(xgd(:), ygd(:), zgd(:));
% [xgd,ygd,zgd] = meshgrid(x,y,z);        [n,m,k] = size(xgd);
[faces, vertices] = freeBoundary(DT);
vert1 = vertices(faces(:,1),:);
vert2 = vertices(faces(:,2),:);
vert3 = vertices(faces(:,3),:);

% vmin = [x(1),y(1),z(1];vmax = [x(2),y(2),z(2)];  
% vert_all = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
% faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
% figure, h = patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
% set(h,'FaceAlpha',0.5);

%% find ray intersect pts with ALL plane faces
% (More than 6 planes of cube, as code uses traingular mesh)
[Intersect,Dist,bcd1,bcd2, Xcord] = TriangleRayIntersection(orig, dir, ...
                                            vert1, vert2, vert3);

%% Find entry and exit to cube, in physical space
% ie find shortest and longest distance of interesects
inter_pos = Xcord(Intersect,:);
inter_dis = Dist(Intersect);

[entry_D,Ien_D] = min(inter_dis);
[exit_D ,Iex_D] = max(inter_dis);

entry_P = inter_pos(Ien_D,:);           % units still in Rs
exit_P  = inter_pos(Iex_D,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Find entry/ exit of cube as index cube value.
% temp = norm(entry_P(1,:));
xn= length(xgdpts); yn = length(ygdpts);
% entry_P = [5.5,3.1,0.2];
ix =  find( (xgdpts-entry_P(1,1))<=0,1,'last') ;
iy =  find( (ygdpts-entry_P(1,2))<=0,1,'last') ;
iz =  find( (zgdpts-entry_P(1,3))<=0,1,'last') ;
Enpt = [ix,iy,iz];
EnIndexGrd = voxel2ind (xn,yn,Enpt);
% ans = [xg(EnIndexGrd),yg(EnIndexGrd),zg(EnIndexGrd)] * (km2rs/1000);

ix =  find( (xgdpts-exit_P(1,1))<=0,1,'last') ;
iy =  find( (ygdpts-exit_P(1,2))<=0,1,'last') ;
iz =  find( (zgdpts-exit_P(1,3))<=0,1,'last') ;
Expt = [ix,iy,iz];
ExIndexGrd = voxel2ind (xn,yn,Expt);


%% Find index array for ray along Data cube
[xpt,ypt,zpt] = bresenham_line3d (Enpt,Expt,0);
brhmPt = [xpt',ypt',zpt'];
BhmIndexGrd = voxel2ind (xn,yn,brhmPt);  % need to check this is same as inbuilt: sub2ind
sub2ind

%% figure checking
if Fplot == 1
    f2 = figure;
    trisurf(faces, vertices(:,1),vertices(:,2),vertices(:,3),Intersect*1.0,'FaceAlpha', 0.7)
    hold on;
    lend = n_f * dir;
    l1 = line('XData',orig(1)+[0 lend(1)],'YData',orig(2)+[0 lend(2)],'ZData',...
      orig(3)+[0 lend(3)],'Color','r','LineWidth',3)
    p1 = plot3(xg(BhmIndexGrd)*(km2rs/1000), yg(BhmIndexGrd)*(km2rs/1000),...
          zg(BhmIndexGrd)*(km2rs/1000), 's','markerface','g')
    axis equal; axis([-10,15, -Inf, Inf, -Inf, Inf])
    xlabel('Xaxis');ylabel('Yaxis');zlabel('Zaxis')
    set(gca, 'CameraPosition', [106.2478  -35.9079  136.4875])
end

%%
OUT_Indices = BhmIndexGrd;
return


