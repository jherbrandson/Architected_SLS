%Control Loop
timesteps = 5;
DefineUnitCells.m;
par = defineUCpar();
 for t = 1:timesteps
%%
% Define Unit Cells of architectured material
% DefineUnitCellCFun(par);
%%
% Simulate Macro problem (beam)
Beamstructure();
%%
% Simulate unit cell Temperature based on IC
SimulateBeamThermal.m;
% Simluate STress on the unit cell based on IC and Thermal Profile
SimulateBeamStress.m;

% Calculate Topological Optimization and Emergy Functions 
%%
% Update Geometry and


% Loop
end