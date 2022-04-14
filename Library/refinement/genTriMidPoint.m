function [Mfine] = genTriMidPoint(Mcoarse)

assert(Mcoarse.n==2,'Mid-point rule is reserved for triangular meshes')

Elt = Mcoarse.elt;
Nelt = size(Elt,1);

[edg2vtx,~,~] = subsimplices(Mcoarse,1);
Nedg = size(edg2vtx,1);
idx = [edg2vtx(:,1); edg2vtx(:,2)];
jdx = [edg2vtx(:,2); edg2vtx(:,1)];
v = [(1:size(edg2vtx,1))';(1:size(edg2vtx,1))'];
VTE = sparse(idx,jdx,v,Nedg,Nedg);

V = Mcoarse.vtx;
Nv = size(V,1);
edgMid = (V(edg2vtx(:,1),:) + V(edg2vtx(:,2),:))/2;
NewV = [V; edgMid]; % Set of vertices of the new mesh.
VTE = VTE+Nv;
% VTE(i,j) = k, where k is the middle of i and j, and i,j, k refer to
% vertices by their indices in V.

ind = 1:Nelt;
NewElem = zeros(Nelt*4,3);
for p = 1:3
    q = mod(p,3) + 1;
    r = mod(q,3) + 1;
    Vp = Elt(:,p);  
    Mpq = VTE((Elt(:,p)-1)*Nedg + Elt(:,q));
    Mpr = VTE((Elt(:,p)-1)*Nedg + Elt(:,r));
    NewElem(4*(ind - 1)+p,:) = [Vp,Mpq,Mpr];
end
M23 = VTE((Elt(:,2)-1)*Nedg + Elt(:,3));
M31 = VTE((Elt(:,3)-1)*Nedg + Elt(:,1));
M12 = VTE((Elt(:,1)-1)*Nedg + Elt(:,2));
NewElem(4*(ind-1)+4,:) = [M23 M31 M12];

neiElt = Mcoarse.nei_elt;
neiFct = Mcoarse.nei_fct;

NewNeiElt = zeros(4*Nelt,3);
NewNeiFct = zeros(4*Nelt,3);

id = [true true true];

for el = 1:Nelt
    for alpha = 1:3
        NewNeiElt(4*(el-1) + 4,alpha) = 4*(el-1) + alpha;
        NewNeiFct(4*(el-1) + 4,alpha) = 1;
        NewNeiElt(4*(el-1) + alpha,1) = 4*(el-1) + 4;
        NewNeiFct(4*(el-1) + alpha,1) = alpha;
    end
    for alpha = 1:3
        qAlpha = mod(alpha,3)+1;
        rAlpha = mod(qAlpha,3)+1;
        Vq = Elt(el,qAlpha);
        Vr = Elt(el,rAlpha);
        
        
        
        nel = neiElt(el,alpha);
        
        if nel ~= 0
            beta = neiFct(el,alpha);
            idbeta = id;
            idbeta(beta) = false;
            qBeta = find(and(idbeta,Elt(nel,:)==Vq));
            rBeta = find(and(idbeta,Elt(nel,:)==Vr));
            
            
            NewNeiElt(4*(el-1) + qAlpha,3) = 4*(nel-1) + qBeta;
            NewNeiElt(4*(el-1) + rAlpha,2) = 4*(nel-1) + rBeta;
            if mod(rBeta-qBeta,3)==1
                NewNeiFct(4*(el-1) + qAlpha,3) = 3;
                NewNeiFct(4*(el-1) + rAlpha,2) = 2;
            else
                NewNeiFct(4*(el-1) + qAlpha,3) = 2;
                NewNeiFct(4*(el-1) + rAlpha,2) = 3;
            end
            
        end
    end
end

Mfine = GeneralizedMesh(NewV,NewElem,NewNeiElt,NewNeiFct);

end

