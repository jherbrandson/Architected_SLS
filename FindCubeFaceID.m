function cubefaceid = FindCubeFaceID(stressmodel)
Totfaces = stressmodel.Geometry.NumFaces;
cubefaceid = 1:6; %initialized
cubefacepts = [0.5 0 0.5;0.5 1 0.5; 0 0.5 0.5; 1 0.5 0.5; 0.5 0.5 0;0.5 0.5 1;];
avpos = zeros(Totfaces,3);
for j = 1: length(cubefaceid)
dist1 = 0.5;
for i = 1:Totfaces %Plot to highlight faces
    C_n = findNodes(stressmodel.Mesh,'Region','Face',i);
    avpos(i,:) = [sum(stressmodel.Mesh.Nodes(1,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(2,C_n))/length(C_n) ...
             sum(stressmodel.Mesh.Nodes(3,C_n))/length(C_n)];
    distT1 = sqrt((cubefacepts(j,1)-avpos(i,1))^2 ...
         +(cubefacepts(j,2)-avpos(i,2))^2 +(cubefacepts(j,3)-avpos(i,3))^2);
    if distT1<dist1
        cubefaceid(j) = i;
        dist1 = distT1;
    end 
    
    %     figure(); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
    %     hold on; scatter3(nodes(1,C_n),nodes(2,C_n),nodes(3,C_n),'r.');  
end
end
end