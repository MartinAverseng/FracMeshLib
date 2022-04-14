function [A] = femAssembly(M,d,Aloc)

[J,Nd] = localToGlobal(M,d);
Nelt = M.nelt;
n = M.n;
Kd = nchoosek(n+1,d+1);
idx = zeros(Kd^2*Nelt,1);
jdx = zeros(Kd^2*Nelt,1);
val = zeros(Kd^2*Nelt,1);

for el = 1:size(J,1)
    
    V = M.vtx(M.elt(el,:),:);
    E = 1:n;
    m_el = msh(V,E);
    Al = Aloc(m_el);
    
    
    
    for k = 1:Kd
        for l = 1:Kd
            ind = k-1 + Kd*(l-1) + Kd^2*(el-1) + 1;
            genFk = J(el,k);
            genFl = J(el,l);
            idx(ind) = genFk;
            jdx(ind) = genFl;
            val(ind) = Al(k,l);
        end
    end
end

A = sparse(idx,jdx,val,Nd,Nd);



end

