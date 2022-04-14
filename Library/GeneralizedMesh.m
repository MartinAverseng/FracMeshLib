classdef GeneralizedMesh
    
    properties
        vtx; % Nvtx x 3 array of reals (coordinates)
        elt; % Nelt x (n+1) array of numbers in {1,...,Nvtx}
        nei_elt; % Nelt x (n+1) array of numbers in {0,..,Nelt}
        nei_fct; % Nelt x (n+1) array of numbers in {0,...,n+1}
    end
    
    methods
        % Class constructor
        function M = GeneralizedMesh(varargin)
            switch nargin
                case 0
                    elt = [];
                    neiElt = [];
                    neiFct = [];
                case 1
                    arg1 = varargin{1};
                    if isa(arg1,'msh')
                        m = arg1;
                        elt = m.elt;
                        vtx = m.vtx;
                        [neiElt,neiFct] = adj(m);
                    end
                case 4
                    vtx = varargin{1};
                    elt = varargin{2};
                    neiElt = varargin{3};
                    neiFct = varargin{4};
                otherwise
                    error('Unavailable case.')
            end
            M.elt = elt;
            M.vtx = vtx;
            M.nei_elt = neiElt;
            M.nei_fct = neiFct;
        end
        
        function nelt = nelt(M)
            nelt = size(M.elt,1);
        end
        function n = n(M)
            n = size(M.elt,2)-1;
        end
        function nvt = nvtx(M)
            nvt = size(M.vtx,1);
        end
        
        function ndvol = ndv(M)
            ndvol = mshNdvolume(M);
        end
        function N = nrm(mesh)
            T = mesh.tgt;
            if (size(mesh,2) == 3)
                N = cross(T{1},T{2},2);
                N = N ./ (sqrt(sum(N.^2,2)) * [1 1 1]);
            elseif (size(mesh,2) == 2) && is2d(mesh)
                N = T * [0 -1 0 ; 1 0 0 ; 0 0 1]';
            else
                error('msh.m : unavailable case')
            end
        end
        function T = tgt(mesh)
            if (size(mesh,2) == 3)
                T = cell(1,3);
                for i = 1:3
                    ip1  = mod(i,3)+1;
                    ip2  = mod(ip1,3)+1;
                    A    = mesh.vtx(mesh.elt(:,ip1),:);
                    B    = mesh.vtx(mesh.elt(:,ip2),:);
                    T{i} = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
                end
            elseif (size(mesh,2) == 2)
                A = mesh.vtx(mesh.elt(:,1),:);
                B = mesh.vtx(mesh.elt(:,2),:);
                T = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
            else
                error('msh.m : unavailable case')
            end
        end
        function s = size(varargin)
            s = size(varargin{1}.elt);
            if (nargin == 2)
                s = s(varargin{2});
            end
        end
        
        % vtx to element matrix
        function vtx2el = vtx2elt(M)
            idx = repelem(1:M.nelt,M.n+1);
            E = (M.elt)';
            jdx = E(:)';
            el2vtx = sparse(idx,jdx,1,M.nelt,M.nvtx);
            vtx2el = el2vtx';
        end
        
        % Subsimplices
        function [F,set2sub,sub2set] = subsimplices(M,d)
            [F,set2sub,sub2set]  = subsetsRows(M.elt,d+1);
        end
        
        function dM = genBoundary(M)
            dM = genBound(M);
        end
        
        % Generalized subfacets
        function [F,gamma,I] = generalizedSubfacets(M,d)
            [F,gamma,I] = genSubFcts(M,d);
        end
        
        % Refine
        function[M] = refine(M,nit)
            switch M.n
                case 1
                    M = genEdgMidPoint(M);
                case 2
                    M = genTriMidPoint(M);
                otherwise
                    error('No tetrahedral refinement implemented')
            end
            if nit > 1
                M = refine(M,nit-1);
            end
        end
        
        
        % Plot
        function plotFracturedMesh(M)
            assert(M.n == 2,'Can only plot triangular fractured meshes');
            % Will only produce good plots for cracked 2D domains.
            m = msh(M.vtx,M.elt);
            plot(m);
            hold on;
            for i=1:M.nelt
                for alpha = 1:M.n+1
                    j = M.nei_elt(i,alpha);
                    if j > 0
                        edgAlpha = setdiff(M.elt(i,:),M.elt(i,alpha));
                        A = M.vtx(edgAlpha(1),:);
                        B = M.vtx(edgAlpha(2),:);
                        C = M.vtx(M.elt(i,alpha),:);
                        O = (A + B + C)/3;
                        I = (A + B)/2;
                        plot3([O(1) I(1)],[O(2) I(2)],[O(3) I(3)],'-green','LineWidth',3);
                        %text(O(1),O(2),O(3),['K_{',num2str(i),'}'],'FontSize',20);
                    end
                end
                
            end
            hold on
            plot3(M.vtx(:,1),M.vtx(:,2),M.vtx(:,3),'bo','MarkerFaceColor','b')
        end
        function dump(M,name,convention)
            
            if ~exist('convention','var')
                convention = 1;
            end
            
            V = M.vtx;
            E = M.elt -1 + convention;
            nE = M.nei_elt -1 + convention;
            nF = M.nei_fct -1 + convention;
            
            fileID = fopen(name,'w');
            fprintf(fileID,'Generalized Mesh\n\n\n');
            fprintf(fileID,'vtx\n');
            printMatrix(V,',','{','}',fileID);
            
            fprintf(fileID,'\n\n');
            fprintf(fileID,'elt\n');
            printMatrix(E,',','{','}',fileID);
            
            fprintf(fileID,'\n\n');
            fprintf(fileID,'neiElt\n');
            printMatrix(nE,',','{','}',fileID);
            
            
            fprintf(fileID,'\n\n');
            fprintf(fileID,'neiFct\n');
            printMatrix(nF,',','{','}',fileID);
            
            fclose(fileID);
        end
    end
    
    
    
    
    
end

