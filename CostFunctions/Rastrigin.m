classdef Rastrigin < Problem
    properties 
        dimensions
    end
    
    methods
        function obj = Rastrigin(d)
            obj.dimensions = d;
        end
        function cost = evaluate(obj,input)
           if (numel(input) ~= obj.dimensions)
               error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
               %NOTE: error seems to destroy object
           else
               A = 10;
               cost = A*obj.dimensions + sum(input.^2) - A.*sum(cos(2*pi.*input));
           end
        end
    end
end
