function [phi] = directAngle(A,B,C)

% Inputs: Three points in R^2
% Returns: the measure in [0,2pi) of the tirgonometric angle phi between 
% the vectors AB and BC around B


zA = A(:,1) + 1i*A(:,2);
zB = B(:,1) + 1i*B(:,2);
zC = C(:,1) + 1i*C(:,2);

phi = mod(angle((zC - zB)./(zA - zB)),2*pi);




end

