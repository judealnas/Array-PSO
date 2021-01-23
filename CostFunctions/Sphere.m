classdef Sphere < Problem
    properties
        dimensions
    end
    
    methods
        function obj = Sphere(d)
            obj.dimensions = d;
        end
        
        function cost = evaluate(obj,input)
            if (numel(input) ~= obj.dimensions)
                error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                %NOTE: error seems to destroy object
            else
                cost = sum(input.^2);
            end
        end
    end
end