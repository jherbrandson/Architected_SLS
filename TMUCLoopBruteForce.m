%Control Loop for Brute force calculation of thermomechanical unit cell sim
% Simulation Settings
appstress = -50e6;  %Applied stress value, in Pa
% HF = 100;%Heat Flux across the Face, W/m2
DeltaT = 500; %
Tinitial = 0;
fixfacept = [0 0.5 0.5];
loadfacept = [1 0.5 0.5];

%% Define Geometry Start
rez1 = 0.05;
%Define Array Centering Points
nudge = 0.225;
c = [0.5, 0, 0.5;nudge, 0, nudge; 1-nudge,0,nudge; 1-nudge,0, 1-nudge;nudge, 0,1-nudge];
ang1 = 0:15:45;
% czero = c(1,:);
% c = c-czero;
% c = c*[cosd(ang1) 0 sind(ang1); 0 1 0; sind(ang1) 0  -cosd(ang1)];%*[cosd(angle1) sind(angle1) 0; sind(angle1) -cosd(angle1) 0; 0 0 1];%
% c = c+czero;
r = [0.09:0.005:0.12];%0.09;%0.1;
theta = 90; 

plotflag = 1;
%%
% [tetconn, tetpoints, bf, p,shp] = DefUCCylArr(c,r(1),theta(1),rez1);
% f = figure(1); f.Color = 'w';
% % trisurf(tetconn,tetpoints(:,1),tetpoints(:,2),tetpoints(:,3))
% plot(shp)
% axis equal

%Define Material Properties
matprop = DefMatProp();

% DefineAllUC
for i = 1:length(r)
    for j = 1: length(ang1)
    
        [tetconn, tetpoints, bf, p,shp] = DefUCCylArr(c,r(i),theta,ang1(j),rez1);

        [thermalresult, thermalmodel]= SimulateUCThermal(tetpoints,tetconn,...
            Tinitial, DeltaT, fixfacept, loadfacept, matprop);
        [stressresult, stressmodel]  = SimulateUCStress(tetpoints,tetconn,...
            appstress, Tinitial, thermalresult, fixfacept, loadfacept, matprop);

        [Emechsum, ETEsum] = CalculateObjective(stressresult,thermalresult,matprop,Tinitial);
        
        Emechs(i,j) = Emechsum;
        ETEs  (i,j) = ETEsum;
        
        if plotflag
            figure(10)
           PlotResults(thermalresult, thermalmodel, stressresult, stressmodel)
        end
    end
end

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