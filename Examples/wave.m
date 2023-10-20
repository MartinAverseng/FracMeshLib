%% Toy script for FEM in space / FD in time applied to the PDE
% Solves 1/c^2 d^2 u/ dt^2 - Delta u = 0 in Omega \ Gamma
% u(0,x) = u0(x)
% du/dt(0,x)= 0
% du/dn = 0 at the boundary of Omega \ Gamma
% Omega regular polygon, Gamma fracture in Omega
% Note that the Neumann condition at the fracture is applied "independently
% on both sides". 

clear all %#ok
close all

cd ..
cd Library
addpath(genpath(pwd));
cd ..
cd Examples

clc;



%% Creating geometry
close all
disp("Creating geometry");

% Domain (disk)
Ndisk = 80;
Omega = mshDisk(Ndisk,1);

% Fracture
V = [Omega.vtx(1,:);Omega.vtx(2,:);Omega.vtx(4,:);Omega.vtx(13,:);Omega.vtx(6,:);Omega.vtx(3,:);Omega.vtx(7,:)];
E = [1 2; 1 3; 3 4; 1 5;3,6;1,7];
mGamma = msh(V,E); 
% Fracture is constructed as a subset of the edges of Omega. 

% Construct the fractured mesh and plot it
M = fracturedMesh(Omega,mGamma);
plotFracturedMesh(M);
hold on;
pGamma = patch('Faces',mGamma.elt,'Vertices',mGamma.vtx,'EdgeColor','r','LineWidth',3,'Marker','.','MarkerSize',25);
axis equal
axis off
title("Generalized mesh $\mathcal{M}^*_{\Omega \setminus \Gamma}$",'Interpreter','latex')

%% Assembling

disp("Assembling")
M = M.refine(2);

ngauss = 7;
domOmega = dom(M,ngauss); % Quadrature rules on elements of M
Lambda0M = GenFem(M,'P1'); % Space of 0-Whitney forms on M 
% (= conforming piecewise linear element in the energy space 
% ||u||^2_{L^2(Omega)} + ||p||^2_{L^2(Omega)} < inf
% where p = weak gradient of u on Omega \ Gamma (not the same as
% distributional gradient)

Mass = integral(domOmega,Lambda0M,Lambda0M);
K = integral(domOmega,grad(Lambda0M),grad(Lambda0M));

%% Time finite-difference scheme and LU factorization


c=1; % wave speed
dt = 0.005; % time-step
theta = 0.26;
delta = 0.53; % Newark (theta,delta) scheme
Op = Mass/dt^2 + c^2*theta*K;

disp('Factorization of operator')
[L,U] = lu(Op);

%% Initial data

posX = 0.2;
posY = 0.65;
rad = 0.2;
U0 = @(X)(((X(:,1)-posX).^2 + (X(:,2) - posY).^2 < rad^2).*(exp(-rad^2./(max(rad^2-((X(:,1)-posX).^2 + (X(:,2) - posY).^2),eps)))));
% bump function supported in B((posX,posY),rad) 
a = integral(domOmega,Lambda0M,U0); % Projection on Finite element space

U0h = Mass\a; 
U1h = U0h; % 0 initial speed

%% Pre-compute the solutions at every time step
disp("Storing solution data")


maxTimeStep = 1000;
pp = cell(maxTimeStep,1);
for n = 1:maxTimeStep
    if (mod(n,50) == 0)
        fprintf('Done up to time step %d\n',n);
    end
    rhs = Mass/dt^2*(2*U1h - U0h) ...
        - c^2*K*(1/2 + delta - 2*theta)*U1h ...
        - c^2*K*(1/2 - delta + theta)*U0h;
    U2h = U\(L\rhs);
    U0h = U1h;
    U1h = U2h;
    pp{n} = U2h;
end

%% Play the movie

close all
figure;
[X,T] = Lambda0M.dof; % Dofs of Lambda0M are given by the generalized vertices. 

disp("Displaying animation")

p = patch('Faces',T,'Vertices',X,'FaceVertexCData',pp{1},'FaceColor','interp','EdgeColor','interp'); 
patch('Faces',mGamma.elt,"Vertices",mGamma.vtx,'LineWidth',3,'EdgeColor','r')

colormap(gray);
caxis([0 0.1])
axis equal;
axis off

for n = 1:maxTimeStep
    p.FaceVertexCData = pp{n};
    drawnow;
end


