%% Define All Unit Cells

rez1 = 0.05;
%% Define Cell 1 - Solid Box
ci = 1;
[x,y,z]=meshgrid([0:rez1:1],[0:rez1:1],[0:rez1:1]);
x = reshape(x, [size(x,1)*size(x,2)*size(x,3),1]);
y = reshape(y, [size(y,1)*size(y,2)*size(y,3),1]);
z = reshape(z, [size(z,1)*size(z,2)*size(z,3),1]);
box = [x,y,z];

shp1 = alphaShape(box);%inf,'HoleThreshold', 0.2);
% [tetconn, tetpoints] =  alphaTriangulation(shp1);
[BF, P] = boundaryFacets(shp1);

ucell{ci}.Points = P;
ucell{ci}.Faces = BF;

%% Define Cell 2 - Cylinder Array
ci = 2;
nudge = 0.2;
hoffsets = [ 0.5, 0, 0.5;nudge 0, nudge;  1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];
angle = 0;
hoffsets = hoffsets*[cosd(angle) 0 sind(angle); sind(angle) 0  -cosd(angle); 0 1 0];%*[cosd(angle) sind(angle) 0; sind(angle) -cosd(angle) 0; 0 0 1];
r = 0.1;
t = 0; 

[tetconn, tetpoints, bf, p] = DefUCCylArr(hoffsets,r,t,rez1);

stlwrite(triangulation(bf, p), 'TestUCS2.stl');
 
ucell{ci}.Points = p;
ucell{ci}.Faces = bf; 

%% Define Cell Type 3
ci = 3;
%Define object array and placement
nudge = 0.2;
hoffsets = [ 0.5, 0, 0.5;nudge 0, nudge;  1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];

[sptsx,sptsy,sptsz] = cylinder(0.1); %sphere(50)
sptsx = reshape(sptsx, [size(sptsx,1)*size(sptsx,2),1]);
sptsy = reshape(sptsy, [size(sptsy,1)*size(sptsy,2),1]);
sptsz = reshape(sptsz, [size(sptsz,1)*size(sptsz,2),1]);
spts = [sptsx sptsy, sptsz];
%Rotate Cylinder
angle = 90;
spts = spts*[cosd(angle) 0 sind(angle); sind(angle) 0  -cosd(angle); 0 1 0];%*[ 0 cosd(angle) sind(angle); 0 sind(angle) -cosd(angle); 1 0 0];%
sconn = delaunayTriangulation(spts);
[K,v] = convexHull(sconn);

Points = [0 0 0];
CL = [1 2 3];
for i  = 1: size(hoffsets,1)
     Points = [Points; sconn.Points+hoffsets(i,:)];
     CL=[CL; K+(i-1)*size(sconn.Points,1)];
end
Points(1,:) = [];
CL(1,:) = [];
sconn2 = triangulation (CL,Points);%sconn1.ConnectivityList,sconn1.Points);

% rez1 = 0.05;
[x,y,z]=meshgrid([0:rez1:1],[0:rez1:1],[0:rez1:1]);
x = reshape(x, [size(x,1)*size(x,2)*size(x,3),1]);
y = reshape(y, [size(y,1)*size(y,2)*size(y,3),1]);
z = reshape(z, [size(z,1)*size(z,2)*size(z,3),1]);
box = [x,y,z];
fv.vertices = Points;
fv.faces = CL;
% fv.faces = fliplr(fv.faces);
while i <=  size(box,1)
    if inpolyhedron(fv,box(i,:),'TOL',0)
        box(i,:) =[];
        i = i-1;
    end
    i = i+1;
end
boxhole= [box;Points];
shp1 = alphaShape(boxhole,0.071);%inf,'HoleThreshold', 0.2);
[tetconn, tetpoints] =  alphaTriangulation(shp1);
[bf, P] = boundaryFacets(shp1);
 stlwrite(triangulation(bf, P), 'TestUCS3.stl');
 
ucell{ci}.Points = P;
ucell{ci}.Faces = bf; 

%% Plot Unit cell
f = figure(1);
f.Color = 'w';
% trisurf(K,sconn.Points(:,1),sconn.Points(:,2),sconn.Points(:,3))
% trisurf(sconn2.ConnectivityList,sconn2.Points)
plot(shp1)
axis equal
% axis off

%%  FUnctions

function [tetconn, tetpoints, bf, p] = DefUCCylArr(c,r,t,rez1)
% nudge = s;
% hoffsets = [ 0.5, 0, 0.5;nudge 0, nudge;  1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];
hoffsets = c;
hoffsets = hoffsets*[cosd(t) 0 sind(t); sind(t) 0  -cosd(t); 0 1 0];
[sptsx,sptsy,sptsz] = cylinder(r); %sphere(50)
sptsx = reshape(sptsx, [size(sptsx,1)*size(sptsx,2),1]);
sptsy = reshape(sptsy, [size(sptsy,1)*size(sptsy,2),1]);
sptsz = reshape(sptsz, [size(sptsz,1)*size(sptsz,2),1]);
spts = [sptsx sptsy, sptsz];
%Rotate Cylinder
angle1 = 90;
spts = spts*[cosd(angle1) 0 sind(angle1); sind(angle1) 0  -cosd(angle1); 0 1 0];
sconn = delaunayTriangulation(spts);
[K,v] = convexHull(sconn);

%Make Cylinder array
Points = [0 0 0];
CL = [1 2 3];
for i  = 1: size(hoffsets,1)
     Points = [Points; sconn.Points+hoffsets(i,:)];
     CL=[CL; K+(i-1)*size(sconn.Points,1)];
end
Points(1,:) = [];
CL(1,:) = [];
sconn2 = triangulation (CL,Points);%sconn1.ConnectivityList,sconn1.Points);

% rez1 = 0.05;
[x,y,z]=meshgrid([0:rez1:1],[0:rez1:1],[0:rez1:1]);
x = reshape(x, [size(x,1)*size(x,2)*size(x,3),1]);
y = reshape(y, [size(y,1)*size(y,2)*size(y,3),1]);
z = reshape(z, [size(z,1)*size(z,2)*size(z,3),1]);
box = [x,y,z];
fv.vertices = Points;
fv.faces = CL;
% fv.faces = fliplr(fv.faces);
while i <=  size(box,1)
    if inpolyhedron(fv,box(i,:),'TOL',0)
        box(i,:) =[];
        i = i-1;
    end
    i = i+1;
end
boxhole= [box;Points];
shp1 = alphaShape(boxhole,0.071);%inf,'HoleThreshold', 0.2);

[tetconn, tetpoints] =  alphaTriangulation(shp1);
[bf, P] = boundaryFacets(shp1);

 stlwrite(triangulation(bf, P), ['UCCyl_',num2str(r(1)),'_',num2str(size(c,1)),'.stl']);
p = P;

end