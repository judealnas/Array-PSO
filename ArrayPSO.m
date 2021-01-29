% Jude Alnas
% University of Alabama

clear
clc
close all

%look in CostFunctions folder in same directory as this .m file for
%implementations of the Problem.m class
addpath ./CostFunctions 
load('patch_patterns.mat')
%%
c = physconst('LightSpeed');
Dk = 2.33; %relative permittivity
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);

nvars = 3; %# of array elements = nvars + 1; 1 assumed fixed at origin
population = 20;
neighborhood_size = 5; %number of particles in a neighborhood

space_lim = 7e-3; %spacing between patch centers; patch dimensions = 3.5mm x 3.5mm; 
global min_position max_position %originally declared global for use in plottng functions
min_position = space_lim; %intrinsically restricts the particle closest to origin 
max_position = 4*lambda;

vel_lims = [lambda/32,lambda/8]; %velocity magnitude min/max
pos_lims = [min_position, max_position]; 

max_it = 500;

%Create cost function
SLLFunc = SLLCost(nvars); %particle positions = element distances from origin
SLLFunc.Dk = Dk;
SLLFunc.freq0 = freq0;
SLLFunc.element_pattern = 10.^(uStripPatch27GHz1x4ArrayPattern.dBGainTotal_1./10);
SLLFunc.theta = uStripPatch27GHz1x4ArrayPattern.Thetadeg;
SLLFunc.search_mode = "valley";
SLLFunc.optim_mode = "fullpattern";


%gbest parameters - globals for use in updateGlobalBest function
global g_best_cost g_best_pos g_best_particle
g_best_cost = Inf; 
g_best_pos = NaN;

%weights
w_vel = 1; %inertia
w_glob = 0.5; %@(x) x*1; 
w_loc = 0.5; %@(x) 1 - x;
w_per = 1;

weights = [w_vel, w_glob, w_loc, w_per];
%% Swarm Init
%initialize each particle; pass in the number of variables and
%position/velocity limits; 
for i = 1:population
    p = Particle(nvars,pos_lims,vel_lims,space_lim);
    %p.evaluate(SLLFunc); %moved evaluate to main loop
    swarm_arr(i) = p; 
end
updateNeighbors(swarm_arr, neighborhood_size); 
%% Optimization
iteration = 1;

while (iteration <= max_it)
    for p = swarm_arr
        p.evaluate(SLLFunc); 
    end
   
    updateLocalBest(swarm_arr);
    updateGlobalBest(swarm_arr);
    
%     if convergenceCheck(swarm_arr)
%         break;
%     end

%     A = swarmLogEntry(swarm_arr);
%     swarm_log(:,:,iteration) = A;
%     T = array2table(A,'VariableNames', {'X','Y','Cost'})
    
    for p = swarm_arr
        p.updateVelocity(weights, g_best_pos);
        p.updatePosition();
        p.enforceLimits();
    end
    updateNeighbors(swarm_arr,neighborhood_size);
    
    iteration = iteration + 1;
end

%% Plotting
 g_best_pos = [0.021096598402695;0.007074494374363;0.014085880612923];

x = SLLFunc.theta;
%note this function returns AF*EF = rad pattern
y = SLLFunc.getFullPattern(g_best_pos);

hfss_x = Opt_7_1_14_1_21_1.Thetadeg;
hfss_y = Opt_7_1_14_1_21_1.dBGainTotalFreq27GHzPhi900000000000002deg;

%Rectangular plot comparing result, AF, and EF
% figure()
% hold on
% plot(x, y, 'DisplayName', 'Full Pattern')
% plot(x, abs(SLLFunc.getAF(g_best_pos)).^2,'DisplayName','AF^2');
% plot(x, SLLFunc.element_pattern,'DisplayName','EF');
% hold off
% legend
% xlabel("Theta (deg)")
% ylabel("Gain")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rectangular plot comparing result, AF, and EF in dB
% figure()
% title("Decibel Patterns")
% hold on
% plot(x, db(y,'power'),'DisplayName','Full Pattern')
% plot(x, db(abs(SLLFunc.getAF(g_best_pos)).^2,'power'),'DisplayName','AF');
% plot(x, db(SLLFunc.element_pattern,'power'),'DisplayName','EF');
% hold off
% legend
% xlabel("Theta (deg)")
% ylabel("Gain (dB)")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Polar plots comparing result, AF, and EF in dB
figure()
title("Decibel Patterns")
polarplot(deg2rad(x),db(y,'power'), ...
          deg2rad(x),db(abs(SLLFunc.getAF(g_best_pos)).^2,'power'), ...
          deg2rad(x),db(SLLFunc.element_pattern,'power'));
legend('FullPattern','AF','EF')
ax = gca;
ax.ThetaZeroLocation = 'top';
rlim([-50 50])

%polar plots comparing calculated result with HFSS simulation
figure()
polarplot(deg2rad(x),db(y,'power'), ...
          deg2rad(hfss_x),hfss_y);
legend('Calculated','HFSS');
ax = gca;
ax.ThetaZeroLocation = 'top';
rlim([-50 50])


% for i = 1:max_it
%     if i == 2 || i == max_it
%         figure() %put first and last swarm states in separate figures
%     end
%     plotCost(SLLFunc);
%     hold on
%     plotSwarmLogEntry(swarm_log(:,:,i));
%     title(sprintf("Iteration %d",i));
%     axis([min_position max_position min_position max_position]);
%     hold on
%     pause(0.1);
% end
%% Utility Functions
function updateNeighbors(particle_arr, hood_size) %To be method of Swarm class
    for p = particle_arr
        d = [particle_arr.position] - p.position; %get displacement vectors
        dmag = sqrt(sum((d.^2),1)); %calculate displacemen magnitudes; 
        [~,ind] = sort(dmag); %sort; smallest element is 0
        p.neighbors = particle_arr(ind(2:2+hood_size - 1)); 
    end
end
function updateGlobalBest(swarm_arr) %To be method of Swarm class
    global g_best_cost g_best_pos g_best_particle
    [g_best_cost, ind_best] = min([swarm_arr.pbest_cost]);
    g_best_particle = swarm_arr(ind_best);
    g_best_pos = g_best_particle.pbest_pos;
    
end
function updateLocalBest(particle_arr)
    for p = particle_arr
        %OPTION: move these routines into Particle class OR vectorize
        p1 = [p p.neighbors];
        [p.lbest_cost, ind_best] = min([p1.cost]);
        p.lbest_pos = p1(ind_best).position;
    end
end

function b = convergenceCheck()
    b = false; %convergence criterion TBD; for now run as many iterations as allowed    
end

function A = swarmLogEntry(swarm_arr)
    P = [swarm_arr.position]'; %rows = particle position vectors
    C = [swarm_arr.cost]; %rows = particle costs
    C = C(:); 
    A = [P C];
end

function plotSwarmLogEntry(entry)
    A = entry;
    X = A(:,1);
    Y = A(:,2);
    Z = A(:,3);
%     S = swarm_log(:,4,i);
    scatter3(X,Y,Z,36,Z,'MarkerFaceColor','white','LineWidth',1);
end

function plotCost(CostObj)
    %plot Standard Test Functions
    global min_position max_position
    x = min_position:0.1:max_position;
    y = x;
    for i = 1:length(x)
        for j = 1:length(y)
            z(i,j) = CostObj.evaluate([x(i) y(j)]);
        end
    end
    surf(x,y,z','EdgeColor','none','FaceAlpha',0.5)
    view(3)
end

