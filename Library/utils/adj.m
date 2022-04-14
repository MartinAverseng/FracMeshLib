function [neiElt,neiFct] = adj(m)

Nelt = size(m.elt,1);
Nvtx = size(m.vtx,1);
n = size(m.elt,2) - 1;
neiElt = zeros(Nelt,n+1);
neiFct = zeros(Nelt,n+1);
E = m.elt';
for alpha = 1:n+1
    fct = setdiff(1:(n+1),alpha);
    idx = repelem((1:Nelt)',n);
    jdx = E(fct,:); jdx = jdx(:);
    Malpha = sparse(idx,jdx,1,Nelt,Nvtx);
    for beta = 1:n+1
        fct = setdiff(1:(n+1),beta);
        idx = repelem((1:Nelt)',n);
        jdx = E(fct,:); jdx = jdx(:);
        Mbeta = sparse(idx,jdx,1,Nelt,Nvtx);
        [I,J,~] = find(Malpha*Mbeta' == n);
        ind = I~=J;
        I = I(ind);
        J = J(ind);
        neiElt(I,alpha) = J;
        neiFct(I,alpha) = beta;
    end
end


end