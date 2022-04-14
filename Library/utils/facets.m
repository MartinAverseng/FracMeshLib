function [fcts,inc] = facets(m)
switch size(m.elt,2)
    case 4
        [fcts,E] = m.fce;
    case 3
        [fcts,E] = m.edg;
    case 2
        [fcts,E] = m.prt;
    case 1
        fcts = [];
end
jdx = repmat(1:size(E,1),1,size(E,2));
idx = E(:);
inc = sparse(idx,jdx,1,size(fcts.elt,1),size(m.elt,1));
end