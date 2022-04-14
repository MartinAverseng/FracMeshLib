function [J,Nd] = localToGlobal(M,d)

% Computes the local to global function J(i,l). 
% M is a generalized mesh, d is the dimension of the d-subfacet. 

Nelt = size(M.elt,1);
n = size(M.elt,2)-1;
Kd = nchoosek(n+1,d+1);
T = zeros(Nelt,Kd);
J = zeros(Nelt,Kd);


[~,elt2sub,~] = subsimplices(M,d);
[~,gamma,I] = generalizedSubfacets(M,d);
Nd = length(gamma);



for genF = 1:length(gamma)
    S = I(genF);
    L = gamma{genF};
    for l = 1:length(L)
        el = L(l);
        ind = (elt2sub(el,:) == S);
        J(el,ind) = genF;
    end
end



end

