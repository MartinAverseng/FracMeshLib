close all;
clear all;

V = zeros(10,3);

sc = 2;
V(1,:) = [0,0,0];
V(2,:) = sc*[1,0,0];
V(3,:) = sc*[cos(pi/3),sin(pi/3),0];
V(4,:) = V(3,:) + V(2,:);
V(5,:) = V(3,:) - V(2,:);
V(6,:) = V(3,:) + V(3,:) - V(2,:);
V(7,:) = V(3,:) + V(3,:) - V(1,:);
V(8,:) = V(1,:) + V(1,:) - V(2,:);
V(9,:) = V(8,:) + V(5,:) - V(1,:);
V(10,:) = V(9,:) + V(5,:) - V(8,:);

E = [[1 2 3];[1 3 5];[1 5 8];[5 8 9];[5,9,10];[5 10 6];[5 6 3];[3 6 7];[3 7 4];[2 3 4]];

mOmega = msh(V,E);
mGamma = msh(V,[3 5]);

M1= fracturedMesh(mOmega,mGamma);
dump(M1,'M1star.txt',0);


%% Generalized subfacet


e = [[1 2 3];[1 3 5];[5 6 3];[3 6 7];[3 7 4];[2 3 4]];
mOmega = msh(V,e);
mGamma = msh(V,[3 1; 3 2;3 7]);

StarS = fracturedMesh(mOmega,mGamma);
dump(StarS,'StarS.txt',0);



