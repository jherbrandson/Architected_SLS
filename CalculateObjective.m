% Calculate Energy from finite element results
function [Emechsum, ETEsum] = CalculateObjective(stressresult, thermalresult,matprop,Tinitial)
CTE = matprop.CTE;
E = matprop.E;
% stressresults.phi = 
sr = stressresult;
tr = thermalresult;
sxx = sr.Stress.sxx;
Emech = zeros(size(sr.Stress.sxx));
ETE = zeros(size(sr.Stress.sxx));
for i = 1:length(sr.Stress.sxx)
    sigma33{i} = [sr.Stress.sxx(i) sr.Stress.sxy(i) sr.Stress.sxz(i);...
                  sr.Stress.sxy(i) sr.Stress.syy(i) sr.Stress.syz(i);...
                  sr.Stress.sxz(i) sr.Stress.syz(i) sr.Stress.szz(i)];
    sigma61{i} = [sr.Stress.sxx(i) sr.Stress.syy(i) sr.Stress.szz(i)...
                  sr.Stress.sxy(i) sr.Stress.syz(i) sr.Stress.sxz(i)];
           
    epsilon33{i} = [sr.Strain.exx(i) sr.Strain.exy(i) sr.Strain.exz(i);...
                    sr.Strain.exy(i) sr.Strain.eyy(i) sr.Strain.eyz(i);...
                    sr.Strain.exz(i) sr.Strain.eyz(i) sr.Strain.ezz(i)];
            
    epsilon61{i} = [sr.Strain.exx(i) sr.Strain.eyy(i) sr.Strain.ezz(i)...
                    sr.Strain.exy(i) sr.Strain.eyz(i) sr.Strain.exz(i)];
    
    Emech(i) = 0.5*sigma61{i}*epsilon61{i}' ;
    ETE(i)  = E*(CTE*(tr.Temperature(i)-Tinitial))^2;
end

Emechsum = sum(Emech);
ETEsum = sum(ETE);

%% Homogenized Elastic tensor
%Calculate the Homogenized compliance tensor