function [V,elt,indV,indE,K] = removeVertices(vtx,elt,L)

[V,I] = setdiff(vtx,L,'rows');
[V,I,J] = intersect(vtx,V,'rows');

ind = ones(size(elt,1),1);
for alpha = 1:size(elt,2)
    ind = and(ind,ismember(elt(:,alpha),I));
end

indE = find(ind);

elt = elt(ind,:);
K = unique(elt(:));
indV = K;
V = vtx(K,:);
K(K) = 1:length(K);

elt(:) = K(elt(:));



end

