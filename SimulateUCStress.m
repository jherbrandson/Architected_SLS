function [stressresult, stressmodel] = SimulateUCStress(tetpoints,tetconn, ... 
    appstress,Tinitial,thermalresult,fixfacept, loadfacept,matprop)

E = matprop.E; % 200e9; % Elastic Modulus, in Pa
nu = matprop.nu; %0.3; % Poisson's Ratio
CTE = matprop.CTE;%16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
density = matprop.density;%7800;  % Mass Density of Material 1 kg/m^3
% appstress = -50e6;  %Applied stress value, in Pa
% deltaT = 500; 
% Tinitial = 0;

fixface = 2;%'F11'; %ID of the fixed bottom face
loadface = 13;%'F8';  %ID d

Plotflag = 1;
%% Create the Model
stressmodel = createpde('structural','static-solid');

%% Define Geometry
% importGeometry(model,stlfilename);
geometryFromMesh(stressmodel,tetpoints',tetconn');
% geometryFromMesh(model,tetpoints2',tetconn2')

%Plot Geometry
if Plotflag
figure(5)
pdegplot(stressmodel,'FaceLabels','on','VertexLabels','on')%% 'EdgeLabels','on'
axis equal
end
%% Generate Mesh
mesh = generateMesh(stressmodel); %Generate mesh
if Plotflag
figure(6)
pdemesh(mesh)
end
%% Find Faces
Totfaces = stressmodel.Geometry.NumFaces;
fixface  = 1;
loadface = Totfaces;
dist1 = 0.5;
dist2 = 0.5;
for i = 1:Totfaces %Plot to highlight faces
    C_n = findNodes(stressmodel.Mesh,'Region','Face',i);
    avpos = [sum(stressmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(3,C_n))/length(C_n)];
    distT1 =sqrt((fixfacept(1)-avpos(1))^2+(fixfacept(2)-avpos(2))^2 +(fixfacept(3)-avpos(3))^2);
    distT2 =sqrt((loadfacept(1)-avpos(1))^2+(loadfacept(2)-avpos(2))^2 +(loadfacept(2)-avpos(2))^2);
    if distT1<dist1
        fixface = i;
        dist1 = distT1;
    end
    if distT2<dist2
        loadface = i;
        dist2 = distT2;
    end  
    %     figure(); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
%     hold on; scatter3(nodes(1,C_n),nodes(2,C_n),nodes(3,C_n),'r.');  
end

%% Define Material Properties 
structuralProperties(stressmodel,'YoungsModulus',E,'PoissonsRatio',nu,...
    'CTE',CTE,'MassDensity',density);

%% Apply Boundary Conditions
% fixface = [5,10];  %These are the end Face we will Fix
% loadface = [3,8];  %These are the face that will loaded
structuralBC(stressmodel,'Face',fixface,'Constraint','fixed') ;
% structuralBC(model,'Face',fixface,'Constraint','roller') 
% structuralBC(model,'Vertex',10,'Constraint','fixed') 

%% Apply Loads
structuralBoundaryLoad(stressmodel,'Face',loadface,'SurfaceTraction',[0, 0, appstress]);
%structuralBodyLoad(model,'GravitationalAcceleration',GAval)  %Set a G load
% Temperature Load
% Tfinal = Tinitial + deltaT;
% structuralBodyLoad(model,'Temperature',Tfinal) %Setting an intial temp scale

structuralBodyLoad(stressmodel,'Temperature',thermalresult); %Setting an intial temp scale
stressmodel.ReferenceTemperature = Tinitial;         %Set a model Ref Temperature
%% Solve the Model
stressresult = solve(stressmodel); % for stationary structural problems and plot

if Plotflag
figure(7)
pdeplot3D(stressmodel,'ColorMapData',stressresult.Displacement.Magnitude, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
end                    