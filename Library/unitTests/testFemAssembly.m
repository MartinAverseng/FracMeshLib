%% Standard assembly
close all;
clear all;

m = mshDisk(10,1);
M = GeneralizedMesh(m);

d = 0;
n = 2;
Qd = nchoosek(n+1,d+1);
Aloc = @(S)(ones(Qd));

A = femAssembly(M,d,Aloc);
spy(A);



%% Gypsilab style assembly

close all;
clear all;

Omega = mshDisk(15,1);


V = [Omega.vtx(1,:);Omega.vtx(2,:);Omega.vtx(4,:);Omega.vtx(13,:);Omega.vtx(6,:);Omega.vtx(3,:);Omega.vtx(7,:)];
E = [1 2; 1 3; 3 4; 1 5;3,6;1,7];
mGamma = msh(V,E);



M = fracturedMesh(Omega,mGamma);
plotFracturedMesh(M);
M = M.refine(1);
hold on;
pGamma = patch('Faces',mGamma.elt,'Vertices',mGamma.vtx,'EdgeColor','r','LineWidth',3,'Marker','.','MarkerSize',25);


Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Mass = integral(Gamma,gfe,gfe);
Stiff = integral(Gamma,grad(gfe),grad(gfe));

[P,D] = eig(full(Mass\Stiff));
[d,I] = sort(diag(D));

%%

U = Stiff\rhs;
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',U,'FaceColor','interp','EdgeColor','interp');
axis equal
colormap(jet)
hold on
patch('Faces',mGamma.elt,"Vertices",mGamma.vtx,'LineWidth',3,'EdgeColor','k');
axis off


%% Choose eig and plot

neig = 10;
[X,T] = dof(gfe);
% X(:,3) = P(:,I(neig));
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',P(:,I(neig)),'FaceColor','interp','EdgeColor','interp');

hold on
patch('Faces',mGamma.elt,"Vertices",mGamma.vtx,'LineWidth',3,'EdgeColor','k');
axis equal
axis off
colormap(jet)
hold on
edg = M.subsimplices(1);
% patch('Faces',edg,'Vertices',M.vtx,'EdgeColor',[0.5 0.5 0.5]);


%% Circle with cracked radius
clear all;
close all;

mOmega = mshDisk(1000,1);
mOmega.vtx = (mOmega.vtx).*(sum(mOmega.vtx.^2,2));
V = mOmega.vtx(mOmega.vtx(:,2)==0,:);

E = [1:(size(V,1)-1); 2:size(V,1)]';
mGamma = msh(V,E);
M = fracturedMesh(mOmega,mGamma);
% plotFracturedMesh(M);
% hold on
% plot(mGamma,'r');


Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Mass = integral(Gamma,gfe,gfe);
Stiff = integral(Gamma,grad(gfe),grad(gfe));

% Eigenvalue problem
[P,D] = eig(full(Mass\Stiff));
[d,I] = sort(diag(D));

neig = 4;
[X,T] = dof(gfe);
X(:,3) = P(:,I(neig));
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');
% hold on
% mGamma.vtx(:,3) = 0.5;
% plot(mGamma,'r');
alpha = 3/2;
dJ = @(x)(1/2*(besselj(alpha-1,x) - besselj(alpha+1,x)));
disp(dJ(sqrt(d(neig))));
axis equal

% Neumann boundary value problem
f = @(X)(X(:,2));
L = integral(Gamma,gfe,f);
C = integral(Gamma,gfe,@(X)(0*X(:,1)+1));


U = (Stiff+C*C')\L;
[X,T] = dof(gfe);
X(:,3) = U;
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');


%% 

clear all
close all;


load('mOmegaJunction');
load('mGammaJunction');
M = fracturedMesh(mOmegaJunction,mGammaJunction);
figure;
plot(mGammaJunction);

Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Stiff = integral(Gamma,grad(gfe),grad(gfe));

% Neumann boundary value problem
f = @(X)(sin(3*atan2(X(:,2),X(:,1))));
L = integral(Gamma,gfe,f);
C = integral(Gamma,gfe,@(X)(0*X(:,1)+1));


U = (Stiff+C*C')\L;
[X,T] = dof(gfe);
X(:,3) = U;
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');

%% Compare generalized meshes and normal meshes

close all;
clear all;


m = mshDisk(100,1);
M = GeneralizedMesh(m);

Gamma = dom(m,3);
gGamma = dom(M,3);
fe = fem(m,'P1');
gfe = GenFem(M,'P1');

Mass = integral(Gamma,fe,fe);
Stiff = integral(Gamma,grad(fe),grad(fe));

gMass = integral(gGamma,gfe,gfe);
gStiff = integral(gGamma,grad(gfe),grad(gfe));

[P1,D1] = eig(full(Mass\Stiff));
[P2,D2] = eig(full(gMass\gStiff));

[d1,I1] = sort(diag(D1));
[d2,I2] = sort(diag(D2));

neig = 4;

m.vtx(:,3) = P1(:,I1(neig));
m.col = P1(:,I1(neig));
patch('Faces',m.elt,'Vertices',m.vtx,'FaceVertexCData',P1(:,I1(neig)),'FaceColor','interp');

[X,T] = dof(gfe);
X(:,3) = P2(:,I2(neig));
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',P2(:,I2(neig)),'FaceColor','interp');


