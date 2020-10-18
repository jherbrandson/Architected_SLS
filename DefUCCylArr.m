function [tetconn, tetpoints, bf, p,shp, phic,phiv] = ...
    DefUCCylArr(c,r,theta,ang1,rez1)
% Define a unit cell topology consisting of an array of cylindrical holes
% c are the centering points, r are the radii (single or array) t is
% rotation angle of cylindrical holes cf,cn,intf,intn,
% nudge = 0.25;
% c = [ 0.5, 0, 0.5;nudge 0, nudge;  1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];
%Rotate Array
czero = c(1,:);
c = c-czero;
c = c*[cosd(ang1) 0 sind(ang1); 0 1 0; sind(ang1) 0  -cosd(ang1)];%*[cosd(angle1) sind(angle1) 0; sind(angle1) -cosd(angle1) 0; 0 0 1];%
c = c+czero;
% Define Cylinder
[sptsx,sptsy,sptsz] = cylinder(r); %sphere(50)
sptsx = reshape(sptsx, [size(sptsx,1)*size(sptsx,2),1]);
sptsy = reshape(sptsy, [size(sptsy,1)*size(sptsy,2),1]);
sptsz = reshape(sptsz, [size(sptsz,1)*size(sptsz,2),1]);
spts = [sptsx sptsy, sptsz];
%Rotate Cylinder
% theta = theta;
spts = spts*[cosd(theta) 0 sind(theta); sind(theta) 0  -cosd(theta); 0 1 0];
sconn = delaunayTriangulation(spts);
[K,v] = convexHull(sconn);

%Make Cylinder array
CPoints = [0 0 0];
CL = [1 2 3];
for i  = 1: size(c,1)
     CPoints = [CPoints; sconn.Points+c(i,:)];
     CL=[CL; K+(i-1)*size(sconn.Points,1)];
end
CPoints(1,:) = [];
CL(1,:) = [];
% sconn2 = triangulation (CL,CPoints);%sconn1.ConnectivityList,sconn1.Points);

% rez1 = 0.05;
[x,y,z] = meshgrid([0:rez1:1],[0:rez1:1],[0:rez1:1]);
x = reshape(x, [size(x,1)*size(x,2)*size(x,3),1]);
y = reshape(y, [size(y,1)*size(y,2)*size(y,3),1]);
z = reshape(z, [size(z,1)*size(z,2)*size(z,3),1]);
box = [x,y,z];
% Test Polyhedron Cylinders
testpolyc.vertices = CPoints;
testpolyc.faces = CL;
% fv.faces = fliplr(fv.faces);


j = 1;
while i <=  size(box,1)
    if inpolyhedron(testpolyc,box(i,:),'TOL',0)
        hole(j,:) =  box(i,:);
        box(i,:) =[];
        i = i-1;
        j = j+1;
    end
    i = i+1;
end
boxhole= [box;CPoints];
phic = [box;CPoints;hole];
phiv = zeros(size(phic,1));
%Assign Phi;
for i  = 1: size(box,1)
    phiv(i) = 1;
end
for i = 1: size(hole,1)
    phiv(end-(i-1)) = -1;
end

shp = alphaShape(boxhole,0.071);%inf,'HoleThreshold', 0.2);

[tetconn, tetpoints] =  alphaTriangulation(shp);
[bf, p] = boundaryFacets(shp);

%%
% %Split into two reference surfaces, the inside and the outside of the cube
% ind1 = zeros(size(CPoints,1));
% insp = CPoints;
% outsp = zeros(size(p,1)-size(CPoints,1),3);
% m = 1; 
% n = 1;
% for i = 1:size(p,1)
%       ind1(i)  = find(all(CPoints' == p(i,:)'),1);
%       if isempty(ind1(i))
%           outsp(n,:) = p(i,:);
%           indp(i) =n;
%           n = n+1;
%            
%       else 
%           insp(m,:) = p(i,:);
%           indp(i) = m;
%           m = m+1;
%       end
% %     indins(i) = find(all(p'==CPoints(i,:)'));
% %     insp(i,:) = p(indins(i),:);
% end
% %Connection of inside and outside surfaces
% m = 1;
% n = 1;
% for j = 1:size(bf,1)
%     if all(bf(j,:)==ind1)
%         for i = 1: size(bf,2)
%             
%             insc(m,i) = indp(bf(j,i));%%
%         end
%         m = m+1;
%     else
%         for i = 1: size(bf,2)
%             outsc(n,i) = indp(bf(j,i));%%
%         end
%         n+1;               
%     end
% end
%%
stlwrite(triangulation(bf, p), ['UCCyl_',num2str(r(1)),'_',num2str(size(c,1)),'.stl']);

end