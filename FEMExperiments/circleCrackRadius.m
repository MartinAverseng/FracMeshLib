%% Circle with cracked radius
clear all;
close all;

N = 500;
mOmega = mshDisk(N,1);
% mOmega.vtx = (mOmega.vtx).*(sum(mOmega.vtx.^2,2));
V = mOmega.vtx(mOmega.vtx(:,2)==0,:);

E = [1:(size(V,1)-1); 2:size(V,1)]';
mGamma = msh(V,E);
M = fracturedMesh(mOmega,mGamma);
% figure;
% plotFracturedMesh(M);
% title('Fractured mesh: circle with cracked radius')

%% Assembly

Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Mass = integral(Gamma,gfe,gfe);
Stiff = integral(Gamma,grad(gfe),grad(gfe));

%% Eigenvalue problem

[P,D] = eig(full(Stiff),full(Mass));
[d,I] = sort(diag(D));

neig = 2;
[X,T] = dof(gfe);
X(:,3) = P(:,I(neig));
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');
axis equal
title('Eigenvalue problem')
view(0,16);

disp(eigenDiskCutRadius());
disp(sqrt(d(1:10)));
%% Neumann boundary value problem
% 
f = @(X)(X(:,2));
L = integral(Gamma,gfe,f);
C = integral(Gamma,gfe,@(X)(0*X(:,1)+1));


U = (Stiff+C*C')\L;
[X,T] = dof(gfe);
X(:,3) = U;
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');
axis equal
title('Neumann BVP')
view(5,0)

%% Neumann BVP with graded mesh
N = 500;
mOmega = mshDiskGraded(N,1);
V = mOmega.vtx(mOmega.vtx(:,2)==0,:);

E = [1:(size(V,1)-1); 2:size(V,1)]';
mGamma = msh(V,E);
M = fracturedMesh(mOmega,mGamma);


Gamma = dom(M,7);
gfe = GenFem(M,'P1');

Mass = integral(Gamma,gfe,gfe);
Stiff = integral(Gamma,grad(gfe),grad(gfe));

L = integral(Gamma,gfe,f);
C = integral(Gamma,gfe,@(X)(0*X(:,1)+1));


U = (Stiff+C*C')\L;
[X,T] = dof(gfe);
X(:,3) = U;
figure
patch('Faces',T,'Vertices',X,'FaceVertexCData',X(:,3),'FaceColor','interp');
axis equal
title('Neumann BVP, graded mesh')
view(5,0)