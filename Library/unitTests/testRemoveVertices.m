close all;
clear all;
m = shuffle(mshSquare(20,[1,1]));
L = m.vtx(1:3,:);
m = shuffle(m);
plot(m);
hold on
plot(L(:,1),L(:,2),'ro');
[vtx,elt,indV,indE] = removeVertices(m.vtx,m.elt,L);
m2 = msh(vtx,elt);
figure;
plot(m2);
close all;
disp('success');