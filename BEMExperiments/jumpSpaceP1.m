function [P,JumpOp] = jumpSpaceP1(M)
% Returns the matrix of the linear jump operator [] : Vh -> Yh,
% where Vh is the discrete multi-trace space on M and Yh a complement of
% the single-trace space.

[F,~,~] = M.generalizedSubfacets(0);

Nvtx = max(F);

Nd = size(F,1);
Ntilde = Nd - Nvtx;
P = zeros(Ntilde,Nd);

k = 0; % number of lines of P filled
for v = 1:Nvtx
    idx = F == v;
    gfs = find(idx);
    for i = 1:(sum(idx)-1)
        P(k+1,gfs(i)) = 1;
        P(k+1,gfs(end)) = -1;
        k = k+1;
    end
end

JumpOp = 0*P;
k = 0; % number of lines of JumpOp filled
for v = 1:Nvtx
    idx = F == v;
    gfs = find(idx);
    L = sum(idx);
    for i = 1:(L-1)
        for j = 1:(L-1)
            if (i==j)
                JumpOp(k+i,gfs(j)) = 1 - 1/L;
            else
                JumpOp(k+i,gfs(j)) = -1/L;
            end
        end
        if (i==L)
            JumpOp(k+i,gfs(L)) = 1/L;
        else
            JumpOp(k+i,gfs(L)) = -1/L;
        end
    end
    k = k+(L-1);
end



end

