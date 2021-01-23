classdef Giddings < Problem
    properties
        dimensions
    end
    
    methods
        function obj = Giddings(d)
            obj.dimensions = d;
        end
        
        function cost = evaluate(obj, input)
            if (numel(input) ~= obj.dimensions)
                    error("Dimension Mismatch: Cost function dimension and particle position dimension do not match");
                    %NOTE: error seems to destroy object
            else
                cost = 0;
                for i = 2:numel(input)
                    cost = cost + sin(i + input(i))*sqrt(abs(cos(input(i-1))*i+1));
                end
            end
        end
    end
end