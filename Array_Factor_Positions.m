%%
c = physconst('LightSpeed');
N = 3; %# of elements
Ai = ones(N,1); %Feeding Coefficient
Dk = 2.33; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant
theta = (-0:0.05:360)'; %degree angle
d = 0.803*lambda;   %uniform element spacing in meters
% d = 5.56e-3
D = (0:1:N-1)*d;   %array of element positions
% D = [0; 0.014277499428709;0.007014985734191;0.029036206480868]
%% Calculate AF
Exp_arr = exp(1i*k*D(:)'.*sin(deg2rad(theta(:))));
AF = Exp_arr*Ai;
AF2 = abs(AF).^2;
%normalize
AF2_max = max(AF2);
ind_max = find(abs(AF2_max - AF2) < 1E-10); 
%% Peak Finding
%Locate indices of peaks; remove index of main lobe; calculate max of remaining peaks

%max SLL using valley finder; valleys used to find main node null points
ind_vall = find(islocalmin(AF2)); %find indices of valleys 
n = sort([ind_vall(:); ind_max(:)]); %sorted array of indices of valleys and the main peaks
null1 = n(find(n == ind_max(1)) + 1); %ind_max(1) = index of boresight
null2 = n(find(n == ind_max(2)) - 1); %ind_max(2) = index of rear boresight

[max_SLL1,ind_max_SLL] = max(AF2(null1:null2))

%max SLL using peak finder
ind_peak = find(islocalmax(AF2)); %find indices of peaks (detection may cause problems)

ind_side = setdiff(ind_peak,ind_max); %remove index of boresights; gives indices of side lobe peaks

[max_SLL2, ind2_max_SLL] = max(AF2(ind_side)) %ind2 = second-order index, an index to retrieve ANOTHER index
max_SLL2 == AF2(ind_side(ind2_max_SLL))

%% Plot
%data to plot located SL peaks
y = zeros(size(theta));
y(ind_side) = AF2_max;

y1 = zeros(size(theta));
y1(ind_peak) = AF2_max;

y2 = zeros(size(theta));
y2(ind_vall) = AF2_max;

figure()
plot(theta, (AF2),theta,y, theta, y1,theta,y2);

figure()
polarplot(deg2rad(theta), abs(AF));
%rlim([-20 0])
