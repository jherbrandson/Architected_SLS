
stlfilename = 'TestBeam1a.stl';
E = 200e9; % Elastic Modulus, in Pa
nu = 0.3; % Poisson's Ratio
CTE1 = 16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
density1 = 7800;  % Mass Density of Material 1 kg/m^3
TC = 16; %Thermal Conductivity, W/mK
Cp = 500; %Specific heat, J/kgK
appstress = -50e6;  %Applied stress value, in Pa
HF = 100;%Heat Flux across the Face, W/m2
DeltaT = 500; %
Tinitial = 0;
Tfinal = Tinitial+DeltaT;
Plotflag = 1;

fixface = 2;%'F11'; %ID of the fixed bottom face
loadface = 13;%'F8';  %ID d

%% Create the Model
thermalmodel = createpde('thermal','steadystate');

%% Define Geometry
% importGeometry(model,stlfilename);
geometryFromMesh(thermalmodel,tetpoints',tetconn')
% geometryFromMesh(model,tetpoints2',tetconn2')

%Plot Geometry
if Plotflag
figure(5)
pdegplot(thermalmodel,'FaceLabels','on','VertexLabels','on')%% 'EdgeLabels','on'
axis equal
end
%% Generate Mesh
mesh = generateMesh(thermalmodel); %Generate mesh
if Plotflag
figure(6)
pdemesh(mesh)
end
%% Define Material Properties 
% structuralProperties(model,'YoungsModulus',E,'PoissonsRatio',nu,...
%     'CTE',CTE1,'MassDensity',density1);
thermalProperties(thermalmodel,'ThermalConductivity',TC,...
    'MassDensity',density1,'SpecificHeat',Cp)
%% Apply Boundary Conditions
% fixface = [5,10];  %These are the end Face we will Fix
% loadface = [3,8];  %These are the face that will loaded
thermalBC(thermalmodel, 'Face',fixface, 'Temperature',Tinitial);
thermalBC(thermalmodel, 'Face',loadface, 'Temperature',Tfinal);

% thermalBC(thermalmodel, 'Face',loadface, 'HeatFlux',HF);

% internalHeatSource(thermalmodel, 'Node', 21)
%% Apply Loads
% thermalIC(thermalmodel, )
% thermalmodel.ReferenceTemperature = Tinitial;         %Set a model Ref Temperature
%% Solve the Model
thermalresult = solve(thermalmodel); % for stationary structural problems and plot
if Plotflag
figure(7)
pdeplot3D(thermalmodel,'ColorMapData',thermalresult.Temperature)%, ...
%                           'Deformation',thermalresult.Displacement, ...
%                           'DeformationScaleFactor',1)

end