% Create a Multilayer Composite Beam Simulation with 2 different materials
% Dr. Ryan Smith, MATE Calpoly cc2019
submatdims = [6 4 5]; %
L = ones([1,submatdims(1)]); % Beam Length, m
W = ones([1,submatdims(2)]);  % Beam Width, m
H = ones([1,submatdims(3)]); % Beam Heights (Thicknesses), m


nu1 = 0.33;  %Poisson's Ratio, Material 1
E1 = 200e9;  %Young's Modulus, Material 1
nu2 = 0.25;  %Poisson's Ratio, Material 2
E2 = 10e6;   %Young's Modulus, Material 2
appstress = [ 1 0 0]; % Applied stress on the end of the beam (z-direction

xcut = 0; % Xvalue of the data section cut
rez = 0.0002;%max(H)/2;%0.002;  %Resolution of the solution cut plot

%%  
%Create Geometry arrays



%%  Create And Solve FEM Structural Static Model of Composite Beam
%Create Model
model = createpde('structural','static-solid');
%Define Geometry
gm = multicuboid(L,W,H,'ZOffset',offsetsz); %Define a rectangular Prism
model.Geometry = gm;
pdegplot(model,'FaceAlpha',0.5,'Facelabels','on')
% geom_T = importGeometry(model_T,'RectPrism.stl'); %Importing geometry to the PDE model object
mesh = generateMesh(model);%,'Hmax', meshrez); %Generate mesh

%% Define Material Properties 
for i = 1:Layers
   if mod(i,2)
structuralProperties(model,'Cell',i,'YoungsModulus',E1,'PoissonsRatio',nu1);
   else 
  structuralProperties(model,'Cell',i,'YoungsModulus',E2,'PoissonsRatio',nu2);     
   end
end
%% Apply Boundary Conditions
fixface = [5,10,15,20,25];%,30];
loadface = [3,8,13,18,23];%,28];
structuralBC(model,'Face',fixface,'Constraint','fixed') 
% structuralBC(model,'Face',fixface,'Constraint','roller') 
% structuralBC(model,'Vertex',10,'Constraint','fixed') 
%% Apply Loads
structuralBoundaryLoad(model,'Face',loadface,'SurfaceTraction',appstress)
% structuralBodyLoad(structuralmodel,'GravitationalAcceleration',GAval)
% structuralBodyLoad(structuralmodel,'Temperature',Tval)

%% Solve the Model
result = solve(model); % for stationary structural problems
%Make a Data cut and Interpolate
yqs = [-W/2:rez:W/2];
if ~mod(Layers,2)
    zqs = [0:rez:sum(H)*(Layers/2)];
else
    zqs = [0:rez:(sum(H)*((Layers-1)/2) + H(1))];
end
[yq zq] = meshgrid(yqs,zqs);
xq = xcut*ones(size(zq));
% uintrp = interpolateSolution(result,xq,yq,zq);
intrpStress = interpolateStress(result,xq,yq,zq);
intrpStrain = interpolateStrain(result,xq,yq,zq);
%% Plot
figure(2)
pdeplot3D(model,'ColorMapData',result.Stress.sxx,'Deformation',result.Displacement)%, 'ElementLabels','on')
figure(3)
sxx = reshape(intrpStress.sxx,size(yq));
px = pcolor(yq,zq,sxx);
px.EdgeColor='none';
caxis1= caxis;
colorbar

%% Stress Plots
figure(4)
sgtitle(['Stress Components at x = ',num2str(xcut)])
pxx = subplot(2,3,1);
% title(['\sigma_x_x at x = ',num2str(xcut)])
posxx = pxx.Position;
sxx = reshape(intrpStress.sxx,size(yq));
px = pcolor(yq,zq,sxx);
px.EdgeColor='none';
caxis1 = caxis;
colorbar
pxx.Position = posxx;
pxx = subplot(2,3,2);
% title(['\sigma_y_y at x = ',num2str(xcut)])
posxx = pxx.Position;
syy = reshape(intrpStress.syy,size(yq));
px = pcolor(yq,zq,syy);
px.EdgeColor='none';
caxis = caxis1;
colorbar
pxx.Position = posxx;
set(gca, 'YTick', []) 
pxx =subplot(2,3,3);
% title(['\sigma_z_z at x = ',num2str(xcut)])
posxx = pxx.Position;
szz = reshape(intrpStress.szz,size(yq));
px = pcolor(yq,zq,szz);
px.EdgeColor='none';
caxis = caxis1;
set(gca, 'YTick', []) 
colorbar
pxx.Position = posxx;
pxx =subplot(2,3,4);
% title(['\sigma_x_y at x = ',num2str(xcut)])
posxx = pxx.Position;
sxy = reshape(intrpStress.sxy,size(yq));
px = pcolor(yq,zq,sxy);
px.EdgeColor='none';
caxis = caxis1;
colorbar
pxx.Position = posxx;
pxx =subplot(2,3,5);
set(gca, 'YTick', []) 
% title(['\sigma_y_z at x = ',num2str(xcut)])
posxx = pxx.Position;
syz = reshape(intrpStress.syz,size(yq));
px = pcolor(yq,zq,syz);
px.EdgeColor='none';
caxis = caxis1;
colorbar
set(gca, 'YTick', []) 
pxx.Position = posxx;
pxx =subplot(2,3,6);
% title(['\sigma_x_z at x = ',num2str(xcut)])
posxx = pxx.Position;
sxz = reshape(intrpStress.sxz,size(yq));
px = pcolor(yq,zq,sxz);
caxis = caxis1;
colorbar
pxx.Position = posxx;
px.EdgeColor='none';
set(gca, 'YTick', []) 
%% Strain Plots
figure(5)
sgtitle(['Strain Components at x = ',num2str(xcut)])
pxx =subplot(2,3,1);
% title(['\epsilon_x_x at x = ',num2str(xcut)])
posxx = pxx.Position;
sxx = reshape(intrpStrain.exx,size(yq));
px = pcolor(yq,zq,sxx);
px.EdgeColor='none';
caxis1 = caxis;
colorbar
pxx.Position = posxx;
pxx = subplot(2,3,2);
% title(['\epsilon_y_y at x = ',num2str(xcut)])
posxx = pxx.Position;
syy = reshape(intrpStrain.eyy,size(yq));
px = pcolor(yq,zq,syy);
px.EdgeColor='none';
caxis = caxis1;
 set(gca, 'YTick', []) 
colorbar
pxx.Position = posxx;
pxx = subplot(2,3,3);
% title(['\sigma_z_z at x = ',num2str(xcut)])
posxx = pxx.Position;
szz = reshape(intrpStrain.ezz,size(yq));
px = pcolor(yq,zq,szz);
px.EdgeColor='none';
caxis = caxis1;
 set(gca, 'YTick', []) 
colorbar
pxx.Position = posxx;
pxx =subplot(2,3,4);
% title(['\epsilon_x_y at x = ',num2str(xcut)])
posxx = pxx.Position;
sxy = reshape(intrpStrain.exy,size(yq));
px = pcolor(yq,zq,sxy);
px.EdgeColor='none';
colorbar
pxx.Position = posxx;
caxis = caxis1;
pxx =subplot(2,3,5);
% title(['\epsilon_y_z at x = ',num2str(xcut)])
posxx = pxx.Position;
syz = reshape(intrpStrain.eyz,size(yq));
px = pcolor(yq,zq,syz);
px.EdgeColor='none';
caxis = caxis1;
set(gca, 'YTick', []) 
colorbar
pxx.Position = posxx;
pxx = subplot(2,3,6);
% title(['\epsilon_x_z at x = ',num2str(xcut)])
posxx = pxx.Position;
sxz = reshape(intrpStrain.exz,size(yq));
px = pcolor(yq,zq,sxz);
caxis = caxis1;
set(gca, 'YTick', []) 
px.EdgeColor='none';
colorbar
pxx.Position = posxx;