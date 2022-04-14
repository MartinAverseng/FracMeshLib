function [MGamma] = extrinsicInflation(Gamma,Omega)


M1 = genBound(fracturedMesh(Omega,Gamma));
dOmega = bnd(Omega);

[vtx,elt,~,indE] = removeVertices(M1.vtx,M1.elt,dOmega.vtx);

K = 0*(1:size(M1.elt,1));
K(indE) = 1:length(indE);
neiElt = K(M1.nei_elt(indE,:));
neiFct = M1.nei_fct(indE,:);

MGamma = GeneralizedMesh(vtx,elt,neiElt,neiFct);


end

