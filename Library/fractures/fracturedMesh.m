function M = fracturedMesh(Omega,Gamma)

[Facets,inc] = facets(Omega);


[Gam,I] = intersect(Facets,Gamma);
assert(size(Gam.elt,1)==size(Gamma.elt,1),'Gamma must be a union of facets of Omega');

M = GeneralizedMesh(Omega);
for i = 1:size(Gamma.elt,1)
    l = find(inc(I(i),:));
    assert(length(l)<=2,'Omega must be a regular mesh');
    if length(l)==2
        alpha = find(M.nei_elt(l(1),:) == l(2));
        beta = M.nei_fct(l(1),alpha);
        M.nei_elt(l(1),alpha) = 0;
        M.nei_fct(l(1),alpha) = 0;
        M.nei_elt(l(2),beta) = 0;
        M.nei_fct(l(2),beta) = 0;
    end
end






end
