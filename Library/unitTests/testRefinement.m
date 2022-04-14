close all;
clear all;

Omega = mshSquare(15,[1 1]);


V = [0 0 0; -0.25,0,0;0.25,0,0;0,0.25,0;0,-0.25,0];
E = [1 2; 1 3; 1 4; 1 5];
mGamma = msh(V,E);

M = fracturedMesh(Omega,mGamma);
tic;
M = refine(M,3);
toc;
plotFracturedMesh(M);
hold on
plot(mGamma,'r');


Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Mass = integral(Gamma,gfe,gfe);
Stiff = integral(Gamma,grad(gfe),grad(gfe));

[P,D] = eig(full(Mass\Stiff));
[d,I] = sort(diag(D));

neig = 9;
[X,T] = dof(gfe);
X(:,3) = P(:,I(neig));
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',P(:,I(neig)),'FaceColor','interp');
hold on
mGamma.vtx(:,3) = 0.5;
plot(mGamma,'r');