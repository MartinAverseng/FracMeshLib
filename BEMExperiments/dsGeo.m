close all;
clear all;

m = ds();
figure
plot(m);
axis equal;
title('Initial mesh of the fracture')

M = intrinsicInflation(m);
M = M.refine(2);
[F,gamma,I] = M.generalizedSubfacets(0);

m = m.midpoint.midpoint;
A = accumarray(F,1);
figure
plot(m.edg,ones(size(m.vtx,1),1)*[0.8 0.8 0.8]);
hold on
ind1 = A==1; 
plot3(M.vtx(ind1,1),M.vtx(ind1,2),M.vtx(ind1,3),'ro','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',7);
ind2 = A==2; 
plot3(M.vtx(ind2,1),M.vtx(ind2,2),M.vtx(ind2,3),'ko','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);

ind3 = A==3;
plot3(M.vtx(ind3,1),M.vtx(ind3,2),M.vtx(ind3,3),'bo','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10);
ind4 = A==4;
plot3(M.vtx(ind4,1),M.vtx(ind4,2),M.vtx(ind4,3),'go','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',12);
ind7 = A==7;
plot3(M.vtx(ind7,1),M.vtx(ind7,2),M.vtx(ind7,3),'yo','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',15);
axis equal;
title('Generalized vertices of the mesh')