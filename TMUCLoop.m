%Control Loop

%% Simulation Settings
appstress = -50e6;  %Applied stress value, in Pa
% HF = 100;%Heat Flux across the Face, W/m2
DeltaT = 500; %
Tinitial = 0;
fixfacept = [0 0.5 0.5];
loadfacept = [1 0.5 0.5];
timesteps = 1; % Total Refinement 
tinc = 0.001;
Vmax= 0.001;

%% Define Geometry Start
rez1 = 0.05; %FEM Simulation Resolution
%Define Array Centering Points
nudge = 0.24;
c = [0.5, 0, 0.5;nudge, 0, nudge; 1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];
ang1 = 0;%0:15:45; %Array Rotation Angle
r = 0.1; %[0.09:0.005:0.12];%0.09;%0.1;  % Radius of the cylinders (passthrough profile?)
theta = 0;  %Cylinder rotation angle (needs centering and proper 

plotflag = 1;
%% Test
[tetconn,tetpoints,bf,p,shp,phic,phiv] = DefUCCylArr(c,r,theta,ang1,rez1);
f = figure(1); f.Color = 'w';
plot(shp) 
axis equal % trisurf(tetconn,tetpoints(:,1),tetpoints(:,2),tetpoints(:,3))

%%
%Define Material Properties
matprop = DefMatProp();

% DefineAllUC
Vt = 0;
phicold = phic;
phicnew = phic;
phivold = phiv;
phivnew = phiv;
for i = 1:length(r)
    for j = 1: length(ang1)
        %Define the starting model, cylinder hole array
        [tetconn,tetpoints,bf,p,shp,phic,phiv] = DefUCCylArr(c,r,theta,ang1,rez1);
        Vt = zeros(size(phic));
        [thermalresult, thermalmodel]= SimulateUCThermal(tetpoints,tetconn,...
            Tinitial, DeltaT, fixfacept, loadfacept, matprop);
        for t = 1: timesteps
        phicold = phicnew;
        
        [stressresult, stressmodel]  = SimulateUCStress(tetpoints,tetconn,...
            appstress, Tinitial, thermalresult, fixfacept, loadfacept, matprop);

        [Emechsum, ETEsum] = CalculateObjective(stressresult,thermalresult,matprop,Tinitial);
        teststress = 1;
%         DH = CalcHomTen(tetpoints,tetconn, phic, teststress, matprop);
        
        Emechs(i,j) = Emechsum;
        ETEs  (i,j) = ETEsum;
        
        f  = Emechsum + ETEsum;
        g = sum(phiv);
        Vt = 
        phicnew = phicold - tinc*Vt;
        
            if plotflag
                figure(10)
               PlotResults(thermalresult, thermalmodel, stressresult, stressmodel)
            end
        
        end
    end
end

%% Find internal Face midpoints
cubefaceid = FindCubeFaceID(stressmodel);
% Totfaces = stressmodel.Geometry.NumFaces;
% cubefaceid = 1:6; %initialized
% cubefacepts = [0.5 0 0.5;0.5 1 0.5; 0 0.5 0.5; 1 0.5 0.5; 0.5 0.5 0;0.5 0.5 1;];
% 
% for j = 1: length(cubefaceid)
% dist1 = 0.5;
% for i = 1:Totfaces %Plot to highlight faces
%     C_n = findNodes(stressmodel.Mesh,'Region','Face',i);
%     avpos(i,:) = [sum(stressmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
%              sum(stressmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
%              sum(stressmodel.Mesh.Nodes(3,C_n))/length(C_n)];
%     distT1 = sqrt((cubefacepts(j,1)-avpos(i,1))^2 ...
%          +(cubefacepts(j,2)-avpos(i,2))^2 +(cubefacepts(j,3)-avpos(i,3))^2);
%     if distT1<dist1
%         cubefaceid(j) = i;
%         dist1 = distT1;
%     end 
%     
%     %     figure(); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
%     %     hold on; scatter3(nodes(1,C_n),nodes(2,C_n),nodes(3,C_n),'r.');  
% end
% end

%% Plot results
if plotflag
    if size(Emechs,1)>1
        if size(Emechs,2)>1
            figure(11)
            pcolor(Emechs)
        elseif size(Emechs,2)==1
            figure(12)
            subplot(2,1,1)
            
                plot(r,Emechs)
                title('Mechanical Energy Sum')  
            subplot(2,1,2)
                plot(r,ETEs)
                title('Thermomech Energy Sum')
        end
    end
end

%% Accesory Functions
function null = PlotResults(thermalresult, thermalmodel, stressresult, model)
%% Plot
figure(1)
pdeplot3D(thermalmodel,'ColorMapData',thermalresult.Temperature)%, ...

figure(2)
pdeplot3D(model,'ColorMapData',stressresult.Displacement.Magnitude, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      
figure(3)
subplot(2,3,1)
pdeplot3D(model,'ColorMapData',stressresult.Stress.sxx, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off
subplot(2,3,2)
pdeplot3D(model,'ColorMapData',stressresult.Stress.syy, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off
subplot(2,3,3)
pdeplot3D(model,'ColorMapData',stressresult.Stress.szz, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off
subplot(2,3,4)
pdeplot3D(model,'ColorMapData',stressresult.Stress.sxy, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off
subplot(2,3,5)
pdeplot3D(model,'ColorMapData',stressresult.Stress.syz, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off
subplot(2,3,6)
pdeplot3D(model,'ColorMapData',stressresult.Stress.sxz, ...
                          'Deformation',stressresult.Displacement, ...
                          'DeformationScaleFactor',1)
                      colorbar off                    
end

function matprop = DefMatProp()
matprop.E = 200e9; % Elastic Modulus, in Pa
matprop.nu = 0.3; % Poisson's Ratio
matprop.CTE = 16e-6;     % Coefficient of Thermal Expansion, Material 1, 1/K
matprop.density = 7800;  % Mass Density of Material 1 kg/m^3
matprop.K = 16; %Thermal Conductivity, W/mK
matprop.Cp = 500; %Specific heat, J/kgK

end