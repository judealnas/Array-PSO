classdef Particle < handle
    properties
        cost = Inf
        dims
        position
        velocity
        
        pbest_cost = Inf
        pbest_pos
        
        lbest_cost = Inf
        lbest_pos
        
        v_lims  %vel_lims = [min max] mag per dimension
        p_lims  %pos_lims = [min max] per dimension
        spacing_lim %minimum spacing between particle dimensions
        
        constraints = []
        
        neighbors Particle
    end
    
    methods 
        function obj = Particle(d,pos_lims, vel_lims, spacing_lim) %later add means of constructing particles from file
            po_lims = sort(pos_lims);
            ve_lims = sort(vel_lims);
            
            %assign object properties
            obj.p_lims = po_lims;
            obj.v_lims = ve_lims;
            obj.dims = d;
            obj.spacing_lim = spacing_lim;
            
            %po_lims = [min pos, max pos]
            obj.position = obj.getRandPos();   
           
            %vel_lims = [minimum vel mag, max vel mag]
            obj.velocity = rand(d,1).*(ve_lims(2)-ve_lims(1))+ vel_lims(1);
            
            obj.enforceLimits();
        end
    
        function evaluate(obj,CostObj)
            obj.cost = CostObj.evaluate(obj.position);
            if (obj.cost < obj.pbest_cost)
                obj.pbest_cost = obj.cost;
                obj.pbest_pos = obj.position;
            end
        end
        
        function updateVelocity(obj, weights, gbest_pos)
            % weights = [wv, wg, wl, wp]
            w = weights(:);
            dg = gbest_pos - obj.position;
            dl = obj.lbest_pos - obj.position;
            dp = obj.pbest_pos - obj.position;
            
            x = [obj.velocity(:) dg(:) dl(:) dp(:)]; %
            r = rand(size(x,2),1); %every weight has random coefficient
            r(1) = 1; %do not randomize velocity
            v = x*(r.*w);
            obj.velocity = v;
        end
        
        function updatePosition(obj)
            obj.position = obj.position + obj.velocity;
        end
        
        function updateLocalBest(obj) %better to vectorize than use this
            p = [obj obj.neighbors];
            [obj.lbest_cost, ind_best] = min([p.cost]);
            obj.lbest_pos = p(ind_best).position;
        end
          
        function enforceLimits(obj) 
            obj.enforcePosLimits();
            obj.enforceVelLimits();
        end
        
        function enforcePosLimits(obj)
            %enforce position limits; currently limit enforced as
            %a ceiling; enforced on individual dimensions 
            lims = sort(obj.p_lims);
            p = obj.position;
            if any(p < lims(1)) || any(p > lims(2)) || any(diff(sort(p)) < obj.spacing_lim)
                p = obj.getRandPos();
            end     
            obj.position = p;
        end    
                
        function enforceVelLimits(obj)
            %%enforce velocity limits; currently limit enforced as
            %%a ceiling; enforced on individual dimensions
            lim = sort(obj.v_lims);
            v = obj.velocity;
            
            %replace element magnitudes above lim(2) with lim(2)
            v(abs(v) > lim(2)) = sign(v(abs(v) > lim(2)))*lim(2);
            
            %replace element magnitudes below lim(1) with lim(1)
            v(abs(v) < lim(1)) = sign(v(abs(v) < lim(1)))*lim(1);
            
            obj.velocity = v;
        end
        
        function pos = getRandPos(obj)
            % Return a randomly-generatged array of positions satisfying position and spacing
            % limits
            p = rand(obj.dims,1).*(max(obj.p_lims) - min(obj.p_lims)) + min(obj.p_lims);
            difference = diff(sort(p)); %get space between consecutive elements
            while (any(abs(difference) < obj.spacing_lim)) 
                p = rand(obj.dims,1).*(max(obj.p_lims) - min(obj.p_lims)) + min(obj.p_lims);
                difference = diff(sort(p));
            end
            pos = p;
        end
        
        function plot(obj, CostObj)
            CostObj.plot(obj.position);
        end
    end
end



