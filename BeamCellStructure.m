%% Beam Cell Structure
DefineAllUC

%Define the substructure matrix for a beam
submatdims = [6 4 5]; %
submat = ones(submatdims);
submat(:,:, 1) = 1;
submat(:,:,5) = 1;
submat(:,:,2) = 2;
submat(:,:,4) = 2;
submat(:,:,3) = 3;

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