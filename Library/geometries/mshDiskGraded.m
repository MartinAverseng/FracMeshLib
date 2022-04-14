function mesh = mshDiskGraded(N,rad)


dr = sqrt(pi*rad^2/N);
dr = (rad/ceil(rad/dr));
r  = (dr:dr:rad);
r = r.^2;

% Angular uniform discretization
rho = cell(length(r),1); theta = rho;
for ir = 1:length(r)
    Ntheta = ceil(sqrt(N)*sqrt(r(ir)));
    tmp = linspace(0,2*pi,Ntheta);
    theta{ir} = tmp(1:end-1)';
    rho{ir}   = r(ir)*ones(length(theta{ir}),1);
end

% Carthesian coordinates
[x,y] = pol2cart(cell2mat(theta),cell2mat(rho));
X     = [0 0 ; x y];

% Unicity test
tmp = unique(X,'rows','stable');
if (max(abs(X-tmp)) > 1e-12)
    error('mshDisk : non unicity of vertices')
end
   
% Delaunay triangulation
DT = delaunayTriangulation(X(:,1),X(:,2));

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(vtx,elt);
end
