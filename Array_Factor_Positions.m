%%
c = physconst('LightSpeed');
N = 4; %# of elements
Ai = ones(N,1); %Feeding Coefficient
Dk = 2.33; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant
theta = (-0:0.05:360)'; %degree angle
d = 1*lambda/2;   %uniform element spacing in meters
d = 5.56e-3
D = (0:1:N-1)*d;   %array of element positions

%% Calculate AF
Exp_arr = exp(1i*k*D(:)'.*sin(deg2rad(theta(:))));
AF = Exp_arr*Ai;

%normalize
AF_max = max(AF);
ind_max = find(abs(AF_max - AF) < 1E-10);
AF = AF./AF_max; 
%% Peak Finding
%Locate indices of peaks; remove index of main lobe; calculate max of remaining peaks

%max SLL using valley finder; valleys used to find main node null points
ind_vall = find(islocalmin(abs(AF))); %find indices of valleys 
n = sort([ind_vall(:); ind_max(:)])
null1 = n(find(n == ind_max(1)) + 1); %ind_max(1) = index of boresight
null2 = n(find(n == ind_max(2)) - 1); %ind_max(2) = index of rear boresight

max_SLL1 = max(AF(null1:null2))

%max SLL using peak finder
ind_peak = find(islocalmax(abs(AF))); %find indices of peaks (detection may cause problems)

ind_side = setdiff(ind_peak,ind_max); %remove index of boresights; gives indices of side lobe peaks

[max_SLL2, ind2_max_SLL] = max(AF(ind_side)) %ind2 = second-order index, an index to retrieve ANOTHER index
max_SLL2 == AF(ind_side(ind2_max_SLL))

%% Plot
%data to plot located SL peaks
y = zeros(size(theta));
y(ind_side) = 1;

y1 = zeros(size(theta));
y1(ind_peak) = 1;

y2 = zeros(size(theta));
y2(ind_vall) = 1;

figure()
plot(theta, (abs(AF)),theta,y, theta, y1,theta,y2);

figure()
polarplot(deg2rad(theta), abs(AF))
%rlim([-20 0])
