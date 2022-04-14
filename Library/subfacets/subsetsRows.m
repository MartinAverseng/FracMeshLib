function [S,set2sub,sub2set] = subsetsRows(T,d)



N = size(T,1);
n = size(T,2);

assert(n<=4,'No conventions for n > 4')
assert(and(1<=d,d<=n),'d must be between 1 and n');


% Conventions on facets order
switch n
    case 1
        switch d
            case 1
                facets = {1};
        end
    case 2
        switch d
            case 1
                facets = {1,2};
            case 2
                facets = {1:2};
        end
    case 3
        switch d
            case 1
                facets = {1,2,3};
            case 2
                facets = {[2,3],[3,1],[1,2]};
            case 3
                facets = {1:3};
        end
    case 4
        switch d
            case 1
                facets = {1,2,3,4};
            case 2
                facets = {[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]};
            case 3
                facets = {[2 3 4],[3 4 1],[4 1 2],[1 2 3]};
            case 4
                facets = {1:4};
        end
end
Qd = length(facets); %Qd = nchoosek(n,d)

Snon_unique = zeros(N*Qd,d);
for alpha = 1:Qd
    ind = (alpha-1)*N + (1:N);
    Snon_unique(ind,:) = T(:,facets{alpha});
end

tmp           = sort(Snon_unique,2);
[S,~,set2sub] = unique(tmp,'rows');

set2sub = reshape(set2sub,N,Qd);

jdx = set2sub(:);
idx = repmat((1:N)',Qd,1);
sub2set = sparse(idx,jdx,1,N,size(S,1));



end

