mOmega = mshSquare(15,[1 1]);
M = GeneralizedMesh(mOmega);

V = [0 0 0; -0.25,0,0;0.25,0,0;0,0.25,0;0,-0.25,0];
E = [1 2; 1 3; 1 4; 1 5];
mGamma = msh(V,E);
M3 = fracturedMesh(mOmega,mGamma);
dM3 = M.genBoundary;


dump(M3,'Mstar3.txt',0);
dump(dM3M3,'dMstar3.txt',0);