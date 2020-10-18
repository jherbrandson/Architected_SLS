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

fv.vertices = spts;
fv.faces = sconn.ConnectivityList;
% fv.faces = fliplr(fv.faces);
testpoint = [0 0 0];
a = inpolyhedron(fv,testpoint,'TOL',0.05)

figure (3)
trisurf(K,sconn.Points(:,1),sconn.Points(:,2),sconn.Points(:,3))
hold on
scatter3(testpoint(1),testpoint(2),testpoint(3))