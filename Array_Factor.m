clear
clc
close all
%%
c = physconst('LightSpeed');
N = 4; %# of element9
Ai = ones(N,1); %Feeding Coefficient
sub_dims = []; %antenna [x y] dimensions
Dk= 1; %substrate relative permittivity/dielectric constant
freq0 = 28e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant
theta = -180:0.05:180; %degree angle
d = 1*lambda/2*ones(1,N);   %uniform element spacing in meters
D = [0 d];   %array of element spacings
%% Calculate AF
AF = 0;
for v = 1:N
    L = sum(D(1:v));
    AF = AF + Ai(v).*exp(1i*k*L*sin(deg2rad(theta)));
end

%normalize
AF_max = max(AF);
ind_max = find(abs(AF_max - AF) < 1E-10);
AF = AF./AF_max; 
%% Peak Finding
%Locate indices of peaks; remove index of main lobe; calculate max of remaining peaks
ind_peak = find(islocalmax(abs(AF))); %find indices of peaks (detection may cause problems)
%ind_vall = find(islocalmin(abs(AF))); %    
ind_side = setdiff(ind_peak,ind_max); %remove index of max peak; gives indices of side lobe peaks

[max_SLL, ind2_max_SLL] = max(AF(ind_side)); %ind2 = second-order index, an index to retrieve ANOTHER index
max_SLL == AF(ind_side(ind2_max_SLL));

%% Plot
%data to plot located SL peaks
y = zeros(size(theta));
y(ind_side) = 1;

y1 = zeros(size(theta));
y1(ind_peak) = 1;

figure()
plot(theta, abs(AF),theta,y, theta, y1);

figure()
polarplot(deg2rad(theta), abs(AF))
%rlim([-20 0])
