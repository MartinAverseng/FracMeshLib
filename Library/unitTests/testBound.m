%% Normal mesh

m = mshCube(10,[1 1 1]);
M = GeneralizedMesh(m);
dM1 = M.genBoundary;
dm = m.bnd;
dM2 = GeneralizedMesh(dm);


%% Fractured mesh


Omega = mshSquare(15,[1 1]);


V = [0 0 0; -0.25,0,0;0.25,0,0;0,0.25,0;0,-0.25,0];
E = [1 2; 1 3; 1 4; 1 5];
Gamma = msh(V,E);



M = fracturedMesh(Omega,Gamma);
plotFracturedMesh(M);
hold on
plot(Gamma,'r');

dM = genBound(M);
[F,gamma,I] = dM.generalizedSubfacets(0);
multiplicities = accumarray(F,1);
assert(isequal(unique(multiplicities),[1 ;4]));
ddM = genBound(dM);
assert(size(dM.elt,1)==24);
assert(isempty(ddM.elt));

close all;
disp("success");