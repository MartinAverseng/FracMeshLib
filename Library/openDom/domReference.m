function [x,w] = domReference(domain)

% Discrete quadrature of reference domain (segment, triangle or
% tetrahedron)


% Particles mesh
if (size(domain.msh.elt,2) == 1)
    error('domReference.m : no quadrature for point meshes')
    
% Edge mesh
elseif (size(domain.msh.elt,2) == 2)
    [x,w] = domGauss_Legendre1D(domain.gss,0,1);
    
% Triangular mesh
elseif (size(domain.msh.elt,2) == 3)
    if (domain.gss == 1)
        x = [1/3  1/3];
        w = 1;
        
    elseif (domain.gss == 3)
        x = [1/6 1/6 
             2/3 1/6
             1/6 2/3];
        w = [1/3 ; 1/3 ; 1/3];
        
    elseif (domain.gss == 7)
        a = (155-sqrt(15))/1200;
        b = (155+sqrt(15))/1200;
        x = [1/3                1/3
            (6-sqrt(15))/21     (6-sqrt(15))/21
            (6-sqrt(15))/21     (9+2*sqrt(15))/21 
            (9+2*sqrt(15))/21   (6-sqrt(15))/21
            (6+sqrt(15))/21     (6+sqrt(15))/21
            (6+sqrt(15))/21     (9-2*sqrt(15))/21
            (9-2*sqrt(15))/21   (6+sqrt(15))/21];
        w = [9/40 ; a ; a ; a ; b ; b ; b];
        
    elseif domain.gss == 12
        A  = 0.063089014491502;
        B  = 0.249286745170910;
        C  = 0.310352451033785;
        D  = 0.053145049844816;
        P1 = 0.025422453185103;
        P2 = 0.058393137863189;
        P3 = 0.041425537809187;
        x(1:domain.gss,1) = [A 1-2*A A     B 1-2*B B     C D 1-C-D 1-C-D C     D];
        x(1:domain.gss,2) = [A A     1-2*A B B     1-2*B D C C     D     1-C-D 1-C-D];
        w = [P1 P1 P1 P2 P2 P2 P3 P3 P3 P3 P3 P3]'*2;
    
    else
        error('domReference.m : number of Gauss points not supported')
    end
    
% Tetrahedron mesh
elseif (size(domain.msh.elt,2) == 4)
    if (domain.gss == 1)
        x = [1/4 1/4 1/4];
        w = 1;
        
    elseif (domain.gss == 4)
        a = (5-sqrt(5))/20;
        b = (5+3*sqrt(5))/20;
        x = [a a a;
            a a b;
            a b a;
            b a a];
        w = [1/4 1/4 1/4 1/4]';
            
    elseif (domain.gss == 5)
        a = 1/4;
        b = 1/6;
        c = 1/2;
        x = [a a a;
            b b b;
            b b c;
            b c b;
            c b b];
        w = [-4/5 9/20 9/20 9/20 9/20]';
            
    elseif (domain.gss == 15)
        a  = 1/4;
        b1 = (7+sqrt(15))/34; 
        b2 = (7-sqrt(15))/34;
        c1 = (13-3*sqrt(15))/34;
        c2 = (13+3*sqrt(15))/34;
        d  = (5-sqrt(15))/20;
        e = (5+sqrt(15))/20;
        x = [a a a;
            b1 b1 b1;
            b1 b1 c1;
            b1 c1 b1;
            c1 b1 b1;
            b2 b2 b2;
            b2 b2 c2;
            b2 c2 b2;
            c2 b2 b2;
            d d e;
            d e d;
            e d d;
            d e e;
            e d e;
            e e d];
        w1 = 6*(2665-14*sqrt(15))/226800;
        w2 = 6*(2665+14*sqrt(15))/226800;
        w = [48/405 w1 w1 w1 w1 w2 w2 w2 w2 30/567 30/567 30/567 30/567 30/567 30/567]';
        
    else
        error('domReference.m : number of Gauss points not supported')
    end
    
% Unknown type
else
    error('domReference.m : unavailable case')
end
end
