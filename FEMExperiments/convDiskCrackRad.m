% Test convergence eig disk cracked radius:

Ns = [10, 40, 160, 640];
h = [];

Eig = eigenDiskCutRadius()';

for i = 1:length(Ns)
    N = Ns(i);
    mOmega = mshDisk(N,1);
    % mOmega.vtx = (mOmega.vtx).*(sum(mOmega.vtx.^2,2));
    V = mOmega.vtx(mOmega.vtx(:,2)==0,:);
    
    E = [1:(size(V,1)-1); 2:size(V,1)]';
    mGamma = msh(V,E);
    M = fracturedMesh(mOmega,mGamma);
    
    %% Assembly
    
    Gamma = dom(M,7);
    gfe = GenFem(M,'P1');
    
    Mass = integral(Gamma,gfe,gfe);
    Stiff = integral(Gamma,grad(gfe),grad(gfe));
    
    %% Eigenvalue problem
    
    [P,D] = eig(full(Stiff),full(Mass));
    [d,I] = sort(diag(D));
    
    l = mOmega.stp;
    h = [h,l(2)];
    col{i} = abs(Eig - sqrt(d(2:6)))./abs(Eig);
end

