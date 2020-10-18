function DH = CalcHomTen(tetpoints,tetconn, phic, teststress, matprop)

E = matprop.E; % 200e9; % Elastic Modulus, in Pa
nu = matprop.nu; %0.3; % Poisson's Ratio
CTE = matprop.CTE;%16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
density = matprop.density;%7800;  % Mass Density of Material 1 kg/m^3
% appstress = -50e6;  %Applied stress value, in Pa
% deltaT = 500; 
% Tinitial = 0;

fixface = 2;%'F11'; %ID of the fixed bottom face
loadface = 13;%'F8';  %ID d

Plotflag = 0;
%% Create the Model of the architectured material
stressmodel = createpde('structural','static-solid');
refmodel  = createpde('structural','static-solid');
%% Define Geometry
% importGeometry(model,stlfilename);
geometryFromMesh(stressmodel,tetpoints',tetconn');
shp1 = alphaShape(phic,0.071);%inf,'HoleThreshold', 0.2);
[tetconn2, tetpoints2] =  alphaTriangulation(shp1);
geometryFromMesh(refmodel,tetpoints2',tetconn2')

%Plot Geometry
if Plotflag
figure(5)
subplot(1,2,1)
pdegplot(stressmodel,'FaceLabels','on','VertexLabels','on')%% 'EdgeLabels','on'
axis equal
subplot(1,2,2)
pdegplot(refmodel,'FaceLabels','on','VertexLabels','on')%% 'EdgeLabels','on'
axis equal
end
%% Generate Mesh
mesh = generateMesh(stressmodel); %Generate mesh
refmesh = generateMesh(refmodel);
if Plotflag
figure(6)
subplot(1,2,1)
pdemesh(mesh)
subplot(1,2,2)
pdemesh(refmesh)
end
%% Find Faces
Totfaces = stressmodel.Geometry.NumFaces;
cubefaceid = 1:6; %initialized
cubefacepts = [0.5 0 0.5;0.5 1 0.5; 0 0.5 0.5; 1 0.5 0.5; 0.5 0.5 0;0.5 0.5 1;];

for j = 1: length(cubefaceid)
dist1 = 0.5;
for i = 1:Totfaces %Plot to highlight faces
    C_n = findNodes(stressmodel.Mesh,'Region','Face',i);
    avpos(i,:) = [sum(stressmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(3,C_n))/length(C_n)];
    distT1 =sqrt((cubefacepts(j,1)-avpos(1))^2+(cubefacepts(j,2)-avpos(2))^2 +(cubefacepts(j,3)-avpos(3))^2);
    if distT1<dist1
        cubefaceid(j) = i;
        dist1 = distT1;
    end
   
    %     figure(); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
%     hold on; scatter3(nodes(1,C_n),nodes(2,C_n),nodes(3,C_n),'r.');  
end
end

Totfaces = stressmodel.Geometry.NumFaces;
cubefaceidref = 1:6; %initialized
cubefacepts = [0.5 0 0.5;0.5 1 0.5; 0 0.5 0.5; 1 0.5 0.5; 0.5 0.5 0;0.5 0.5 1;];

for j = 1: length(cubefaceidref)
dist1 = 0.5;
for i = 1:Totfaces %Plot to highlight faces
    C_n = findNodes(refmodel.Mesh,'Region','Face',i);
    avpos(i,:) = [sum(refmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
             sum(refmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
             sum(refmodel.Mesh.Nodes(3,C_n))/length(C_n)];
    distT1 =sqrt((cubefacepts(j,1)-avpos(1))^2+(cubefacepts(j,2)-avpos(2))^2 +(cubefacepts(j,3)-avpos(3))^2);
    if distT1<dist1
        cubefaceidref(j) = i;
        dist1 = distT1;
    end
   
    %     figure(); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
%     hold on; scatter3(nodes(1,C_n),nodes(2,C_n),nodes(3,C_n),'r.');  
end
end
%% Define Material Properties 
structuralProperties(stressmodel,'YoungsModulus',E,'PoissonsRatio',nu,...
    'CTE',CTE,'MassDensity',density);
structuralProperties(refmodel,'YoungsModulus',E,'PoissonsRatio',nu,...
    'CTE',CTE,'MassDensity',density);

for i = 1:3 
    for j = 1:3
        for k  = 1:3
            for l  = 1:3
                testload = [0 0 0];
                testload(i) = teststress;
                fixface = cubefaceid(2*j-1);
                fixfaceref = cubefaceid(2*j-1);
                loadface = cubefaceid(2*j);
                loadfaceref = cubefaceid(2*j);
%% Apply Boundary Conditions
structuralBC(stressmodel,'Face',fixface,'Constraint','fixed') ;
structuralBC(refmodel,'Face',fixfaceref,'Constraint','fixed') ;
% structuralBC(model,'Face',fixface,'Constraint','roller') 
% structuralBC(model,'Vertex',10,'Constraint','fixed') 

%% Apply Loads
structuralBoundaryLoad(stressmodel,'Face',loadface,'SurfaceTraction',testload);
structuralBoundaryLoad(refmodel,'Face',loadfaceref,'SurfaceTraction',testload);

%% Solve the Model
stressresult = solve(stressmodel); % for stationary structural problems and plot
refresult = solve(refmodel);

DH(i,j,k,l) =  stressresult(
            end
        end
    end
end 

