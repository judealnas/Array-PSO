%% Encapsulates the weighting coefficients for continuous-domain PSO
classdef Weights < handle 
    
    %Function handles that point to dynamic weighting functions
    properties
        iteration
        w_velocity_func = 1  
        w_glob_func = 1
        w_loc_func = 1
        w_personal_func = 1
        rand_lims = [0 1]
    end
    
    %Coefficients dependent on previous functions
    properties (Dependent)
        w_velocity
        w_glob
        w_loc
        w_personal
    end
           
    methods
        function obj = Weights()
            
        end
        
        function updateWeights(obj,iteration)
            obj.iteration = iteration;
            obj.w_velocity = obj.w_velocity_func(iteration);
            obj.w_glob = obj.w_glob_func(iteration);
            obj.w_velocity = obj.w_loc_func(iteration);
            obj.w_velocity = obj.w_personal_func(iteration);
        end
        
%%%%%%% velocity weight get/seet methods
        function w = get.w_velocity(obj)
            w = obj.randWeight(obj.w_velocity);
        end
        
        function set.w_velocity_func(obj,val)
            obj.w_velocity_func = obj.checkWeightInput(val);
        end

%%%%%%% global best weight get/set methods
        function w = get.w_glob(obj)
            w = obj.randWeight(obj.w_glob);
        end
        
        function set.w_glob_func(obj,val)
            obj.w_glob_func = obj.checkWeightInput(val);
        end

%%%%%%% Local best weight get/set methods
        function w = get.w_loc(obj)
            w = obj.randWeight(obj.w_loc);
        end
        
        function set.w_loc_func(obj,val)
            obj.w_loc_func = obj.checkWeightInput(val);
        end
        
%%%%%%%% Personal best weight get/set methods
        function w = get.w_personal(obj)
            w = obj.randWeight(obj.w_personal);
        end
        
        function set.w_personal_func(obj,val)
            obj.w_personal_func = obj.checkWeightInput(val);
        end
  
    end
    
    methods (Access = private)
        %private function called by get methods to apply random weighting
        function wn = randWeight(obj,w)
            rng('shuffle'); %shuffle seed to approach true randomness
            lims = sort(obj.rand_lims);
            r = lims(1) + (lims(2) - lims(1))*rand(1);
            wn = w*r;
        end
        
        %private function called by set methods to check for constant
        %weight values. if constant create appropriate function handle
        function fh = checkWeightInput(val)
            if isnumeric(val)
                fh = @(i) val;
            elseif isa(val, 'function_handle')
                fh = val;                
            end
        end
    end
end