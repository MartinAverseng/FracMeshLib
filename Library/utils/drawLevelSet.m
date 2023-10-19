function [] = drawLevelSet(M,u,ls,col)

% Reserved for 2d plotting
% M: msh
% u: values at generalized vertices
% l: level
% col: color
assert(M.n==2);
J = localToGlobal(M,0);
locVal = zeros(M.n+1,1);
Xloc = zeros(3,3);
Xbary = zeros(2,3);
for el = 1:M.nelt
    for j = 1:3
        Xloc(j,:) = M.vtx(M.elt(el,j),:);
        locVal(j) = u(J(el,j));
    end
    [sortedVals,I] = sort(locVal);
    for idx = 1:length(ls)
        l = ls(idx);
        if (l>sortedVals(1) && l<sortedVals(end))
            leftBool = l<sortedVals(2);
            for k = 1:2
                if leftBool
                    k1 = 1;
                    k2 = k+1;
                else
                    k1 = k;
                    k2 = 3;
                end
                a = sortedVals(k1);
                b = sortedVals(k2);
                lambda = (l-b)/(a-b);
                Xbary(k,:) = lambda*Xloc(I(k1),:)+(1-lambda)*Xloc(I(k2),:);
            end
            plot([Xbary(1,1),Xbary(2,1)],[Xbary(1,2),Xbary(2,2)],'-','Color',col,'LineWidth',0.5);
        end
    end
end



end

