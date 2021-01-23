classdef Rosenbrock < Problem
    properties
        dimensions
    end
    
    methods
        function obj = Rosenbrock(d)
            obj.dimensions = d;
        end
        function cost = evaluate(obj,input)
            if (numel(input) ~= obj.dimensions)
                error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                %NOTE: error seems to destroy object
            else
                cost = 0;
                for i = 1:numel(input)-1
                    cost = cost + 100*(input(i+1)-(input(i))^2)^2+(1-input(i))^2;
                end
            end
        end
    end
end