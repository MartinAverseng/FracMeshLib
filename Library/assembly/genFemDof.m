function[X,elt2dof] = genFemDof(gfe)

M = gfe.msh; % generalized mesh
n = M.n; % dim of M
V = M.vtx;

switch(gfe.typ)
    case 'P0'
        d = n;
    case 'P1'
        d = 0;
    case 'NED'
        d = 1;
    case 'RWG'
        d = n-1;
end

[F,~,~] = generalizedSubfacets(M,d);
X = 0;
for alpha=1:d+1
    X = X + 1/(d+1)*V(F(:,alpha),:);
end

elt2dof = localToGlobal(M,d);


end