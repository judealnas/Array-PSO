classdef CrossTray < Problem
    properties
        dimensions = 2
    end
    
    methods
        function obj = CrossTray(d)
        end
        
        function cost = evaluate(obj,input)
            if (numel(input) ~= obj.dimensions)
                error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                %NOTE: error seems to destroy object
            else
                x = input(1);
                y = input(2);
                cost = -0.0001*(abs(sin(x)*sin(y)*exp(abs(100-sqrt(x^2+y^2)/pi)))+1)^0.1;
            end
        end
    end
end