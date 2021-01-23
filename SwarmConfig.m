%% Encapsulate paramters of the swarm optimization
classdef SwarmConfig < handle
    properties
        config_name = strcat("Optimization ", string(datetime));
        dimensions
        pos_lims
        vel_lims
        rand_lims = [0 1]
        population
        max_it 
        weights Weights
    end
    
    methods (Static)
        function obj = SwarmConfig(obj, dimensions, pos_lims, vel_lims, population, max_it, weights)
            if (nargin > 0)
                obj.dimensions = dimensions;
                obj.pos_lims = pos_lims;
                obj.vel_lims = vel_lims;
                obj.population = population;
                obj.max_it = max_it;
                obj.weights = weights;
            end
        end
        
        function updateWeights(obj,iteration)
            obj.iteration = iteration;
            obj.weights.updateWeights(iteration);
        end
        
    end
    
end
