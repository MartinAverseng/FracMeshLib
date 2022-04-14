function [dM] = genBound(M)

KF = 0*M.elt;
nb = 0; 

nBound = sum(sum(M.nei_elt == 0));
tmpEl = zeros(nBound,M.n);
nei_elt = zeros(nBound,M.n);
nei_fct = zeros(nBound,M.n);

for i1 = 1:M.nelt
    for alpha1 = 1:M.n+1
        if M.nei_elt(i1,alpha1)==0 
            F1 = M.elt(i1,M.elt(i1,:)~=M.elt(i1,alpha1));
            if KF(i1,alpha1) == 0
                nb = nb+1;
                KF(i1,alpha1) = nb;
                tmpEl(nb,:) = F1;
            end
            for beta1 = 1:M.n
                S = F1(F1 ~= F1(beta1));
                [i2,alpha2] = endChain(M,i1,alpha1,S);
                F2 = M.elt(i2,M.elt(i2,:)~=M.elt(i2,alpha2));
                if KF(i2,alpha2) == 0
                    nb = nb + 1;
                    KF(i2,alpha2) = nb;
                    tmpEl(nb,:) = F2;
                end
                beta2 = find(~ismember(F2,S));
                nei_elt(KF(i1,alpha1),beta1) = KF(i2,alpha2);
                nei_fct(KF(i1,alpha1),beta1) = beta2;
            end
        end
    end
end


V = M.vtx;
I = tmpEl(:);
[Iun,~,K] = unique(I);

elt = tmpEl;
elt(:) = K;
vtx = V(Iun,:);

dM = GeneralizedMesh(vtx,elt,nei_elt,nei_fct);


end

