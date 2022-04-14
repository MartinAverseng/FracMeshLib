function[mT] = mshTorus(N,r1,r2)

N1 = fix(sqrt(r1/r2*N))+1;
N2 = fix(N/N1)+1;

a = N1/N2;
mS = mshSquare(N,[a,1]);


phi = 2*pi/a*mS.vtx(:,1);
theta = 2*pi*mS.vtx(:,2);
r = r1 + r2*sin(theta);
z = r2*cos(theta);

mS.vtx = [r.*cos(phi),r.*sin(phi),z];

mS = swap(mS); % normal must go out.
mT = mshCleanTol(mS);


    function[mesh] =  mshCleanTol(mesh)
        [~,I,J] = uniquetol(single(mesh.vtx),'ByRows',true);
        
        
        % Update mesh
        mesh.vtx = mesh.vtx(I,:);
        if (size(mesh.elt,1) == 1)
            J = J';
        end
        mesh.elt = J(mesh.elt);
        
        % Extract vertex table from element
        Ivtx              = zeros(size(mesh.vtx,1),1);
        Ivtx(mesh.elt(:)) = 1;
        mesh.vtx          = mesh.vtx(logical(Ivtx),:);
        
        % Reorder elements
        Ivtx(Ivtx==1) = 1:sum(Ivtx,1);
        if size(mesh.elt,1) == 1
            mesh.elt = Ivtx(mesh.elt)';
        else
            mesh.elt = Ivtx(mesh.elt);
        end
        
        % Empty mesh
       
            I = (mesh.elt(:,1)==mesh.elt(:,2)) + ...
                (mesh.elt(:,1)==mesh.elt(:,3)) + ...
                (mesh.elt(:,2)==mesh.elt(:,3)) ;
            
            % Tetrahedral mesh
       
        
        % Extract non degenerated elements
        if sum(I)
            mesh = mesh.sub(~I);
        end
    end


end

