function [MGamma] = intrinsicInflation(Gamma)



V = Gamma.vtx;
E = sort(Gamma.elt,2);
n = size(E,2);
nel = size(Gamma.elt,1);
S = [E;fliplr(E)]; % duplicate each element,
% and give them opposite orientation
neiElt = zeros(2*nel,n);
neiFct = zeros(2*nel,n);

for i = 1:size(Gamma.elt,1)
    for alpha = 1:n
        % Get list of simplices that share the facet alpha of
        % Ki
        facet = setdiff(E(i,:),E(i,alpha));
        VJ = [];
        J = [];
        betaJ = [];
        for beta = 1:n
            ind = setdiff(1:n,beta);
            bool = ismember(E(:,ind),facet,'rows');
            bool(i) = 0;
            VJ = [VJ;V(E(bool,beta),:)]; %#ok
            J = [J;find(bool)]; %#ok
            betaJ = [betaJ;bool(bool==1)*beta]; %#ok
        end
        if size(VJ,1)==0
            neiElt(i,alpha) = i + nel;
            neiElt(i+nel,n+1-alpha) = i;
            neiFct(i,alpha) = n+1 - alpha;
            neiFct(i+nel,n+1 - alpha) = alpha;
        else
            % Compute angles between simplices
            
            if mod(alpha,2) == 1
                if n==2
                    a= 2*pi-getAngles(V(facet,:),V(E(i,alpha),:),VJ);
                elseif n==3
                    a = getAngles(V(facet,:),V(E(i,alpha),:),VJ);
                end
            else
                if n==2
                    a = getAngles(V(facet,:),V(E(i,alpha),:),VJ);
                elseif n==3
                    e1 = V(facet(2),:);
                    e2 = V(facet(1),:);
                    P = V(E(i,alpha),:);
                    a = angleTri(e1,e2,P,VJ);
                end
            end
            [~,i1] = min(a); 
            [~,i2] = max(a);
            J1 = J(i1); % The neighbor of i is J1 or J1 + nel
            beta1 = betaJ(i1); % facet is the facet number beta1 of J1
            J2 = J(i2); % The neighbor of i + nel is J2 or J2+nel
            beta2 = betaJ(i2); % facet is the facet number beta2 of J2
            % The one to pick is the one that has an opposite orientation of
            % the facet
            if mod(alpha - beta1,2)==1 % Then opposite orientations of facets
                neiElt(i,alpha) = J1;
                neiElt(J1,beta1) = i;
                neiFct(i,alpha) = beta1;
                neiFct(J1,beta1) = alpha;
            else
                neiElt(i,alpha)= J1 + nel;
                neiElt(J1 + nel,n+1 - beta1) =i;
                neiFct(i,alpha) = n+1 - beta1;
                neiFct(J1 + nel,n+1 - beta1) = alpha;
            end
            if mod(alpha - beta2,2)==1
                % Then facet is oriently oppositely for J1+nel and J2+nel
                neiElt(i+nel,n+1 - alpha)= J2+nel;
                neiElt(J2 + nel,n+1 - beta2) = i+nel;
                neiFct(i+nel,n+1 - alpha) = n+1 - beta2;
                neiFct(J2 + nel,n+1 - beta2) = n+1 - alpha;
            else
                neiElt(i+nel,n+1 - alpha) = J2;
                neiElt(J2,beta2) = i+nel;
                neiFct(i+nel,n+1 - alpha) = beta2;
                neiFct(J2,beta2) = n+1 - alpha;
            end
        end
    end
    
end

MGamma = GeneralizedMesh(V,S,neiElt,neiFct);


end

