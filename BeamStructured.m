
%Define the substructure matrix for a beam
submatdims = [6 4 5]; %
submat = ones(submatdims);
submat(:,:, 1) = 1;
submat(:,:,5) = 1;
submat(:,:,2) = 2;
submat(:,:,4) = 2;
submat(:,:,3) = 3;

rez1 = 0.05;
%% Define Cell 1
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

%% Define Cell 2
ci = 2;
nudge = 0.2;
hoffsets = [ 0.5, 0, 0.5;nudge 0, nudge;  1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];

[sptsx,sptsy,sptsz] = cylinder(0.1); %sphere(50)
sptsx = reshape(sptsx, [size(sptsx,1)*size(sptsx,2),1]);
sptsy = reshape(sptsy, [size(sptsy,1)*size(sptsy,2),1]);
sptsz = reshape(sptsz, [size(sptsz,1)*size(sptsz,2),1]);
spts = [sptsx sptsy, sptsz];
%Rotate Cylinder
angle = 90;
spts = spts*[cosd(angle) 0 sind(angle); sind(angle) 0  -cosd(angle); 0 1 0];
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
% [tetconn, tetpoints] =  alphaTriangulation(shp1);
[bf, P] = boundaryFacets(shp1);
 stlwrite(triangulation(bf, P), 'TestUCS2.stl');
 
ucell{ci}.Points = P;
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

%%
%% Make Full Beam

% %submat = reshape(submat, [size(submat,1)*size(submat,2)*size(submat,3),1]);
allpoints = [0 0 0];
allfaces = [ 1 2 3 ];
m = 1 ;
for i = 1: size(submat,1)
    for j = 1:size(submat,2)
        for k = 1: size(submat,3)
            allfaces  = [allfaces; size(allpoints,1) + ucell{submat(i,j,k)}.Faces ];
            allpoints = [allpoints ; ucell{submat(i,j,k)}.Points+[(i-1),(j-1),(k-1)] ];
           
            cells{m}.Origin = [i, j, k];
            cells{m}.Points = ucell{submat(i,j,k)}.Points+[(i-1),(j-1),(k-1)]; 
            cells{m}.Faces = ucell{submat(i,j,k)}.Faces;
            m = m+1;
        end
    end
end
allpoints(1,:) = [];
allfaces(1,:) = [];
allfaces = allfaces-1;
alltri = triangulation(allfaces,allpoints);
stlwrite(alltri,'TestBeam1a.stl');

shp2 = alphaShape(allpoints,0.078,'HoleThreshold', 4*rez1);%);%0.71; %inf,
[tetconn2, tetpoints2] =  alphaTriangulation(shp2);
[bf2, P2] = boundaryFacets(shp2);
stlwrite(triangulation(bf2, P2), 'TestBeam1b.stl');

%% Plot
f2 = figure(2);
f2.Color = 'w';
trisurf(allfaces,allpoints(:,1),allpoints(:,2),allpoints(:,3))

f3 = figure(3);
f3.Color = 'w';
trisurf(bf2,P2(:,1),P2(:,2),P2(:,3))
% plot(shp2)
axis off