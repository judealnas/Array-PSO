classdef HolderTable < Problem
    properties 
        dimensions = 2
    end
    
    methods
        function obj = HolderTable()
            obj.dimensions = 2;
        end
        function cost = evaluate(obj,input)
           if (numel(input) ~= obj.dimensions)
               error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
               %NOTE: error seems to destroy object
           else
               x = input(1);
               y = input(2);
               cost = -abs(sin(x)*cos(y)*exp(abs(sqrt(x^2+y^2)/pi)));
           end
        end
    end
end
