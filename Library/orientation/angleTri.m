function [a] = angleTri(e1,e2,P,PJ)

% This function takes as arguments an oriented edge e, specified by its two
% endpoints e1 and e2 (the orientation is e1->e2), a point P and a list of
% points P1,...PJ. 

% Construct the triangles T and T1,...,TJ by joining e to P, PJ repsectively. 
% On T, consider the normal vector proportional to the wedge product 
% (e1,e2) x (e1,P). Then, for each triangle Tj, the data e,P,P1 permits to
% define a partition of R^3 into two cylindrical sectors, one of which
% contains the vector n. Then this function compute the angle a(j) of the
% latter.

% We project everything in the plane passing through e1 and perpendicular
% to (e1,e2)
% We then compute the angle in [0,pi] and take 2pi - a when needed. 

J = size(PJ,1);
tau = e2 - e1;
tau = tau/scal(tau,tau);



eP = P - e1;
n = cross(tau,eP);

pieP = (eP - scal(eP,tau)*tau)';
ePj = PJ - ones(J,1)*e1;
s = scal(ePj,n);
piePj = (ePj - scal(ePj,tau)*tau)';
a = zeros(J,1);
for j = 1:J
    a(j) = atan2(norm(cross(pieP,piePj(:,j))), dot(pieP,piePj(:,j)));
end

a(s < 0) = 2*pi - a(s<0);
a(abs(s) < 1e-8) = pi;

function[s] = scal(a,b)

s = a(:,1).*b(:,1) + a(:,2).*b(:,2) + a(:,3).*b(:,3);

end

end
