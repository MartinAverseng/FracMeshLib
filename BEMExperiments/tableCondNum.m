[m,M] = ds();

I = 1:3;
GenCondW = 0*I;
condWTilde = 0*I;
dimKer = 0*I;
condWtilde = 0*I;
for i = I
    F = M.generalizedSubfacets(0);
    s = m.stp;
    W = bemAssembly(M);
    P = jumpSpaceP1(M);
    Wtilde = P*W*P';
    
    
    [~,D] = svd(W);
    [~,Dtilde] = svd(Wtilde);
    
    dimJumps = size(F,1) - max(F);
    dimKer(i) = max(F);
    d = sort(diag(D),'descend');
    dtilde = sort(diag(Dtilde),'descend');
    GenCondW(i) = d(1)/d(dimJumps);
    condWtilde(i) = dtilde(1)/dtilde(end);
    
   
    m = m.midpoint;
    if i<I(end)
        M = M.refine(1);
    end
end


