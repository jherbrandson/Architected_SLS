function [thermalresult, thermalmodel] = SimulateUCThermal(tetpoints,tetconn,...
    Tinitial,DeltaT,fixfacept, loadfacept,matprop)

E = matprop.E; % 200e9; % Elastic Modulus, in Pa
nu = matprop.nu; %0.3; % Poisson's Ratio
CTE = matprop.CTE;%16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
density = matprop.density;%7800;  % Mass Density of Material 1 kg/m^3
K = matprop.K;%16; %Thermal Conductivity, W/mK
Cp = matprop.Cp;%500; %Specific heat, J/kgK
% appstress = -50e6;  %Applied stress value, in Pa
% HF = 100;%Heat Flux across the Face, W/m2
% DeltaT = 500; %
% Tinitial = 0;
Tfinal = Tinitial+DeltaT;

Plotflag = 1;

%fixface = 2;%'F11'; %ID of the fixed bottom face  %loadface = 13;%'F8';  %ID d

%% Create the Model
thermalmodel = createpde('thermal','steadystate');

%% Define Geometry
% importGeometry(model,stlfilename);
geometryFromMesh(thermalmodel,tetpoints',tetconn');
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
%% Find Faces
Totfaces = thermalmodel.Geometry.NumFaces;
fixface  = 1;
loadface = Totfaces;
dist1 = 0.5;
dist2 = 0.5;
for i = 1:Totfaces %Plot to highlight faces
    C_n = findNodes(thermalmodel.Mesh,'Region','Face',i);
    avpos = [sum(thermalmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
             sum(thermalmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
             sum(thermalmodel.Mesh.Nodes(3,C_n))/length(C_n)];
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
% structuralProperties(model,'YoungsModulus',E,'PoissonsRatio',nu,...
%     'CTE',CTE1,'MassDensity',density1);
thermalProperties(thermalmodel,'ThermalConductivity',K,...
    'MassDensity',density,'SpecificHeat',Cp);
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