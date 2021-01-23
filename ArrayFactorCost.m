classdef ArrayFactorCost < handle
    properties
        freq0 %operating frequency
        Dk %relative permittivity
        elements %number of elements
        Ai %feeding coefficients
        dtheta = 0.01
        theta = -180:dtheta:180
    end
    
    methods
        function obj = ArrayFactorCost(el)
            obj.elements = el;
        end
            
        function c = evaluate(obj, positions)
            if (numel(positions) ~= obj.elements)
                c = inf;
            else
                
            end
        end
    end
    
end