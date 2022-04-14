[m,M] = ds();

I = 1:6;
multiTrace = 0*I;
jumps = 0*I;
hs = 0*I;

for i = I
    F = M.generalizedSubfacets(0);
    s = m.stp;
    hs(i) = s(2);
    multiTrace(i) = size(F,1);
    jumps(i) = multiTrace(i) - max(F);
    m = m.midpoint;
    if i < I(end)
        M = M.refine(1);
    end
end

