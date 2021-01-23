classdef McCormick < Problem
    properties
        dimensions = 2
    end
    
    methods
        function obj = McCormick(d)
            
        end
        
        function cost = evaluate(obj, input)
            if (numel(input) ~= obj.dimensions)
                    error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                    %NOTE: error seems to destroy object
            else
                x = input(1);
                y = input(2);
                cost = sin(x+y)+(x-y)^2 - 1.5*x + 2.5*y + 1;
            end
        end
    end
end