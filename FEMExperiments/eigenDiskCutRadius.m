function[z] = eigenDiskCutRadius()

z = [];
check = 0;
for n = 0:10
    alpha = n/2;
    dJ = @(x)(1/2*(besselj(alpha-1,x) - besselj(alpha+1,x)));
    zNew = AllZeros(dJ,1,6,10000);
    if ~isempty(zNew)
        check = max(check,max(abs(dJ(zNew))));
    end
    z = [z,zNew];
end


z = sort(z);
z = z(1:5);



function z=AllZeros(f,xmin,xmax,N)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros
if (nargin<4)
    N=100;
end
options=optimset('Display','off');
options.tolX = 1e-12;
options.tolFun = 1e-12;
dx=(xmax-xmin)/N;
x2=xmin;
y2=f(x2);
z=[];
for i=1:N
    x1=x2;
    y1=y2;
    x2=xmin+i*dx;
    y2=f(x2);
    if (y1*y2<=0)                              % Rolle's theorem : one zeros (or more) present
        z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2),options)]; % Linear approximation to guess the initial value in the [x1,x2] range.
    end
end

end

end