classdef Swarm < handle
    properties
        config %type SwarmConfig object
        
        particle_arr = []
        cost_func
        iteration = 0
        gbest_cost
        gbest_pos
    end
    
    methods
        function Swarm(obj, conf)
        %call other intialization methods here
        end
        
        function defineCostFunc(obj, fhandle)
            obj.cost_func = fhandle;
        end
        
        function configSwarm(obj, pop, max_it, weights, dimensions)
            obj.population = pop;
            obj.max_it = max_it;
            obj.weigts = weights;
            obj.dimensions = dimensions;
        end
        
        function initSwarm(obj)
            import Particle
            for i = 1:obj.population
                particle = Particle();
                obj.particle_arr(i) = particle;
                
            end
            
        end
        function optimize(obj)
            convergence = 0;
            obj.iteration = 0;
            while (1) %do-while structure
                obj.updateWeights();
                
                %determine costs at current particle positions
                for p = obj.particle_arr
                    p.evaluate(obj.cost_func)
                end
                
                %get new global bests
                obj.updateGlobalBests();
                
                %move particles 
                for p = obj.particle_arr
                    p.updateVelocity();
                    p.updatePosition();
                end
                
                obj.updateNeighbors();
                
                convergence = obj.checkConvergence();
                
                if (convergence || (obj.iteration >= obj.config.max_it))
                    break
                end
                
                obj.iteration = obj.iteration + 1;
            end
                        
        end
        
        function c = checkConvergence(obj)
            c = 1;
            for p = obj.particle_arr
                if (p.pbest_position ~= obj.gbest_pos)
                    c = 0;
                end
            end
        end
        
        function updateNeighbors(obj, num_neighbors)
            for p = obj.particle_arr
                d = obj.particle_arr.position - p.position; %get displacement vectors
                dmag = sqrt(sum((d.^2),1)); %calculate displacemen magnitudes
                [~,ind] = sort(dmag); %sort 
                
                %ignore first element as it will always be 0 (vector -
                %itself = 0)
                p.neighbors = obj.particle_arr(ind(2:2+num_neighbors - 1));
            end
        end
        
        function updateGlobalBests(obj)
            old_best_cost = obj.gbest_cost;
            new_best_cost = old_best_cost;
            
            for p = obj.particle_arr
                if (p.cost > new_best_cost)
                    new_best_cost = p.cost;
                    new_best_pos = p.position;
                end
                obj.gbest_cost = new_best_cost;
                obj.gbest_pos = new_best_pos;
            end
        end
        
        function updateWeights(obj)
            obj.config.updateWeights(obj.iteration);
        end
        
    end
end