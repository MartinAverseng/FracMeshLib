%% Test speed 

close all;
clear all;

m = mshCube(1e3,[1 1 1]);
M = GeneralizedMesh(m);
[Fs,labels,idx] = M.generalizedSubfacets(0);


%% Test multi screens

close all;
clear all;

load('mGamma');
MGamma = intrinsicInflation(mGamma);
[F,gamma,I] = generalizedSubfacets(MGamma,0);
disp('success');