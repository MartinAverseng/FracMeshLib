
[m,M] = ds();

I = 1:3;
e1 = 0*I;
e2 = 0*I;
e3 = 0*I;
n1 = 0*I;
n2 = 0*I;
for i = I
    F = M.generalizedSubfacets(0);
    s = m.stp;
    t = tic;
    W = bemAssembly(M);
    toc(t);
    [P,Jump] = jumpSpaceP1(M);
    Wtilde = P*W*P';
    
    
    Gamma = dom(M,7);
    Vh = GenFem(M,'P1'); 
    uinc{1} = @(X)(0*X(:,1));
    uinc{2} = uinc{1};
    uinc{3} = @(X)(0*X(:,1) + 1);
    
    L = integral(Gamma,ntimes(Vh),uinc);
    Ltilde = P*L;
    [x,~,~,it1] = pcg(W,L,1e-6,size(W,1));
    [xtilde,~,~,it2] = pcg(Wtilde,Ltilde,1e-6,size(Wtilde,1));
    n1(i) = it1;
    n2(i) = it2;
    
    e1(i) = norm(W*x - L,2);
    e2(i) = norm(W*P'*xtilde - L,2);
    e3(i) = norm(Wtilde*Jump*x - Ltilde,2);
    if i < I(end)
        M = M.refine(1);
    end
end


