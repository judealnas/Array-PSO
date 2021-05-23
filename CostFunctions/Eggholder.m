classdef Eggholder < Problem
    properties 
        dimensions
    end
    
    methods
        function obj = Eggholder()
            obj.dimensions = 2;
        end
        function cost = evaluate(obj,input)
           if (numel(input) ~= obj.dimensions)
               error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
               %NOTE: error seems to destroy object
           else
               x = input(1);
               y = input(2);
                   
               term1 = -(y+47) * sin(sqrt(abs(y+x/2+47)));
               term2 = -x * sin(sqrt(abs(x-(y+47))));
               
               cost = term1 + term2;
           end
        end
    end
end
