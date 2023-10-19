%% Toy script for FEM at a random fracture network
% Solves -Delta u + 0.01*u = sin(6*pi*Y) in Omega \ Gamma
% du/dn = 0 at the boundary of Omega \ Gamma
% Omega rectangle of side 1, Gamma random fracture in Omega
% Note that the Neumann condition at the fracture is applied "independently
% on both sides". 

clear all %#ok
close all


cd Library
addpath(genpath(pwd));
cd ..

clc;


%% Creating geometry
close all
disp("Creating geometry");

% Domain (disk)
Nmesh = 70;
Omega = mshSquare(Nmesh,[1,1]);

% Fracture
I = rand(size(Omega.edg.elt,1),1)<0.35;
mGamma = Omega.edg.sub(I);
mGamma = setdiff(mGamma,Omega.bnd);
ind = find(ismember(mGamma.vtx,Omega.bnd.vtx,'rows'));
ind = find(sum(ismember(mGamma.elt,ind),2)==0);
mGamma = mGamma.sub(ind);




% Fracture is constructed as a subset of the edges of Omega. 


M = fracturedMesh(Omega,mGamma);
plotFracturedMesh(M);
hold on;
pGamma = patch('Faces',mGamma.elt,'Vertices',mGamma.vtx,'EdgeColor','r','LineWidth',3,'Marker','.','MarkerSize',25);
axis equal
axis off
title("Generalized mesh $\mathcal{M}^*_{\Omega \setminus \Gamma}$",'Interpreter','latex')

%% Assembling

disp("Assembling")
M = M.refine(3);

domOmega = dom(M,7); % Quadrature rules on elements of M
Vh = GenFem(M,'P1'); % Space of 0-Whitney forms on M 
% (= conforming piecewise linear element in the energy space 
% ||u||^2_{L^2(Omega)} + ||p||^2_{L^2(Omega)} < inf
% where p = weak gradient of u on Omega \ Gamma (not the same as
% distributional gradient)


[X,T] = Vh.dof; % Dofs of Vh are given by the generalized vertices. 
Mass = integral(domOmega,Vh,Vh);
K = integral(domOmega,grad(Vh),grad(Vh));

%% 

rhs = integral(domOmega,Vh,@(X)(sin(6*pi*X(:,2))));

%% Solving

U = (Mass + 0.01*K)\rhs;

%% Plotting solution

close all;
figure;
patch('Faces',T,'Vertices',X,'FaceVertexCData',U,'FaceColor','interp','EdgeColor','interp'); 
patch('Faces',mGamma.elt,"Vertices",mGamma.vtx,'LineWidth',3,'EdgeColor','k')
caxis([min(min(U)),max(max(U))]);
colormap(jet);
axis equal
title("Solution u")
hold on

c = colormap;
[c1,c2] = caxis;
nl = size(colormap,1);
l = c1 + (c2 - c1)*(1:nl)/nl;
drawLevelSet(M,U,l,'k');
axis off
title("Plot of solution U")
colorbar

