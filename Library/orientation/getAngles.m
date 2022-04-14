function[a] = getAngles(facet,Vi,VJ)

n = size(facet,1)+1;
if n ==2
    J = size(VJ,1);
    A = ones(J,1)*Vi;
    B = ones(J,1)*facet;
    C = VJ;
    a = directAngle(A,B,C);
elseif n==3
    A = facet(1,:);
    B = facet(2,:);
    a = angleTri(A,B,Vi,VJ);
else
    error("Unavailable case")
end

end