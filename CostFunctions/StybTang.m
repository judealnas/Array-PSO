classdef StybTang < Problem
    properties
        dimensions
    end
    
    methods
        function obj = StybTang(d)
            obj.dimensions = d;
        end
        
        function cost = evaluate(obj,input)
            if (numel(input) ~= obj.dimensions)
                error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                %NOTE: error seems to destroy object
            else
                cost = (sum(input.^4) - 16*sum(input.^2) + 5*sum(input))/2;
            end
        end
    end
end