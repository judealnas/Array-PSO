classdef Griewank < Problem
    properties 
        dimensions
    end
    
    methods
        function obj = Griewank(d)
            obj.dimensions = d;
        end
        function cost = evaluate(obj,input)
           if (numel(input) ~= obj.dimensions)
               error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
               %NOTE: error seems to destroy object
           else
               x = sum(input.^2)/4000;
               y = 1;
               for i = 1:obj.dimensions
                y = y * cos(input(i)/sqrt(i));
               end
               cost = x - y + 1;
           end
        end
    end
end
