%% Simulate Stress on UC

stlfilename = 'TestBeam1a.stl';
E = 200e9; % Elastic Modulus, in Pa
nu = 0.3; % Poisson's Ratio
CTE1 = 16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
density1 = 7800;  % Mass Density of Material 1 kg/m^3
appstress = -50e6;  %Applied stress value, in Pa
deltaT = 500; 
Tinitial = 0;

fixface = 2;%'F11'; %ID of the fixed bottom face
loadface = 13;%'F8';  %ID d

Plotflag = 1;
%% Create the Model
model = createpde('structural','static-solid');

%% Define Geometry
% importGeometry(model,stlfilename);
geometryFromMesh(model,tetpoints',tetconn')
% geometryFromMesh(model,tetpoints2',tetconn2')

%Plot Geometry
if Plotflag
figure(5)
pdegplot(model,'FaceLabels','on','VertexLabels','on')%% 'EdgeLabels','on'
axis equal
end
%% Generate Mesh
mesh = generateMesh(model); %Generate mesh
if Plotflag
figure(6)
pdemesh(mesh)
end
%% Define Material Properties 
structuralProperties(model,'YoungsModulus',E,'PoissonsRatio',nu,...
    'CTE',CTE1,'MassDensity',density1);

%% Apply Boundary Conditions
% fixface = [5,10];  %These are the end Face we will Fix
% loadface = [3,8];  %These are the face that will loaded
structuralBC(model,'Face',fixface,'Constraint','fixed') 
% structuralBC(model,'Face',fixface,'Constraint','roller') 
% structuralBC(model,'Vertex',10,'Constraint','fixed') 

%% Apply Loads
structuralBoundaryLoad(model,'Face',loadface,'SurfaceTraction',[0, 0, appstress])
%structuralBodyLoad(model,'GravitationalAcceleration',GAval)  %Set a G load
% Temperature Load
Tfinal = Tinitial + deltaT;
% structuralBodyLoad(model,'Temperature',Tfinal) %Setting an intial temp scale

structuralBodyLoad(model,'Temperature',thermalresult) %Setting an intial temp scale
model.ReferenceTemperature = Tinitial;         %Set a model Ref Temperature
%% Solve the Model
stressresult = solve(model); % for stationary structural problems and plot
if Plotflag
figure(7)
pdeplot3D(model,'ColorMapData',stressresult.Displacement.Magnitude, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
end                    