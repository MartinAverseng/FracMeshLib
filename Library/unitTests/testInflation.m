%% 2D
clear all; close all;

Omega = mshSquare(15,[1 1]);


V = [0 0 0; -0.25,0,0;0.25,0,0;0,0.25,0;0,-0.25,0];
E = [1 2; 1 3; 1 4; 1 5];
Gamma = msh(V,E);


MGamma1 = extrinsicInflation(Gamma,Omega);
MGamma2 = intrinsicInflation(Gamma);

%% 3D

clear all; close all;
load('mGamma');
load('mOmega');

MGamma1 = extrinsicInflation(mGamma,mOmega);
MGamma2 = intrinsicInflation(mGamma);

[F1,gamma1,I1] = generalizedSubfacets(MGamma1,0);
[F2,gamma2,I2] = generalizedSubfacets(MGamma2,0);

A1 = accumarray(F1,1);
A2 = accumarray(F2,1);


% subplot(1,2,1);
% plot(mGamma.edg,'k');
% hold on
% hold on
% ind1 = A1==1; 
% plot3(MGamma1.vtx(ind1,1),MGamma1.vtx(ind1,2),MGamma1.vtx(ind1,3),'ro','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
% ind3 = A1==3;
% plot3(MGamma1.vtx(ind3,1),MGamma1.vtx(ind3,2),MGamma1.vtx(ind3,3),'bo','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10);
% ind4 = A1==4;
% plot3(MGamma1.vtx(ind4,1),MGamma1.vtx(ind4,2),MGamma1.vtx(ind4,3),'go','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10);
% ind7 = A1==7;
% plot3(MGamma1.vtx(ind7,1),MGamma1.vtx(ind7,2),MGamma1.vtx(ind7,3),'yo','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',20);
% title('extrinsic inflation')
figure
plot(mGamma.edg,'k');
hold on
ind1 = A2==1; 
plot3(MGamma2.vtx(ind1,1),MGamma2.vtx(ind1,2),MGamma2.vtx(ind1,3),'ro','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
ind3 = A2==3;
plot3(MGamma2.vtx(ind3,1),MGamma2.vtx(ind3,2),MGamma2.vtx(ind3,3),'bo','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10);
ind4 = A2==4;
plot3(MGamma2.vtx(ind4,1),MGamma2.vtx(ind4,2),MGamma2.vtx(ind4,3),'go','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',10);
ind7 = A2==7;
plot3(MGamma2.vtx(ind7,1),MGamma2.vtx(ind7,2),MGamma2.vtx(ind7,3),'yo','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',20);
title('intrinsic inflation')
disp('success');
