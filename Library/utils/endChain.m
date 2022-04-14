function [i,tau] = endChain(M,i,alpha,S)

b = true;
while b
    K = M.elt(i,:);
    tau = find(~ismember(K,[K(alpha),S]));
    j = M.nei_elt(i,tau);
    beta = M.nei_fct(i,tau);
    if j == 0
        b = false;
    else
        i = j;
        alpha = beta;
    end
end




end

