t = tic; 
m = mshDisk(100,1);
Gamma = dom(m,7);
Vh = dirichlet(fem(m,'P1'),m.bnd);

G = @(X,Y)femGreenKernel(X,Y,'[1/r]',0);
N = integral(Gamma,Gamma,nxgrad(Vh),G,nxgrad(Vh))/(4*pi);
Nreg = regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh))/(4*pi);
N = N+Nreg;

uinc{1} = @(X)(0*X(:,1));
uinc{2} = uinc{1};
uinc{3} = @(X)(0*X(:,1) + 1);

L = integral(Gamma,ntimes(Vh),uinc);
lambda1 = N\L;


M = intrinsicInflation(m);
GammaInfl = dom(M,7);
MultiTraceSpace = GenFem(M,'P1');

W = bemAssembly(M);
Jump = jumpSpaceP1(M);
Wtilde = Jump*W*Jump';


uinc{1} = @(X)(0*X(:,1));
uinc{2} = uinc{1};
uinc{3} = @(X)(0*X(:,1) + 1);

L = integral(GammaInfl,ntimes(MultiTraceSpace),uinc);
Ltilde = Jump*L;

lambda2 = Jump'*(Wtilde\Ltilde);


% Compare the two: sign differences are due to arbitrary conventions in
% choosing one element out of the two that occupy the same location
% 1e-3 errors are due to sloppy treatment of singular integrals in Gypsilab
disp([lambda1, Jump*lambda2]);
toc(t);