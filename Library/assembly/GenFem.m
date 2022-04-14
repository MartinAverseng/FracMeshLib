classdef GenFem < fem
    
    properties
    end
    
    methods
        function [X,elt2dof] = dof(gfe)
            [X,elt2dof] = genFemDof(gfe);
        end
        
    end
end

