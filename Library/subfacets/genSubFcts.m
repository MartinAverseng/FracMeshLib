function[F,gamma,I] = genSubFcts(M,d)

[subs,~,sub2elt] = subsimplices(M,d);
Nf = size(subs,1);
I = [];
F = [];
ngf = 0;
for i = 1:Nf
    ei = sparse(i,1,1,Nf,1);
    l = find(sub2elt*ei);
    elt_l = M.elt(l,:);
    fct_l = 0*elt_l + 1;
    for j = 1:size(subs,2)
        fct_l(elt_l == subs(i,j)) = 0;
    end
    visited = 0*l;
    visited(:) = false;
    for k = 1:length(l)
        nel = 0;
        if ~visited(k)
            ngf = ngf + 1;
            I(ngf) = i; %#ok
            F(ngf,:) = subs(i,:); %#ok
            gamma{ngf} = []; %#ok
            aux(k);
        end
    end
end

    function[] = aux(k)
        nel = nel+1;
        gamma{ngf}(nel) = l(k);
        visited(k) = true;
        alphas = find(fct_l(k,:));
        for alpha = 1:length(alphas)
            el = M.nei_elt(l(k),alphas(alpha));
            if el ~= 0
                if ~visited(l==el)
                    m = find(l == el);
                    aux(m);
                end
            end
        end
    end
end