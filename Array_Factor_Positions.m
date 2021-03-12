close all

%%
c = physconst('LightSpeed');
N = 6; %# of elements
Ai = ones(N,1); %Feeding Coefficient
Dk = 1; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant
theta = (0:0.04:360)'; %degree angle
d = 1*lambda/2;   %uniform element spacing in meters
% d = 5.56e-3
D = (0:1:N-1)*d;   %array of element positions
% D = [0; 0.014277499428709;0.007014985734191;0.029036206480868]
%% Calculate AF
Exp_arr = exp(1i*k*sin(deg2rad(theta))*D(:)');
AF = Exp_arr*Ai;
AF2 = abs(AF).^2;
AF2(AF2 == 0) = NaN;
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

[max_SLL2, ind2_max_SLL] = max(AF2(ind_side)); %ind2 = second-order index, an index to retrieve ANOTHER index
max_SLL2 == AF2(ind_side(ind2_max_SLL));

%% Plot
%data to plot located SL peaks
y = zeros(size(theta));
y(ind_side) = 1;

y1 = zeros(size(theta));
y1(ind_peak) = 1;

y2 = zeros(size(theta));
y2(ind_vall) = 1;

%Replacing nulls with NaN to eliminate spikes
rlim_low_db = -20;
AF2_norm = AF2/AF2_max;

figure()
plot(theta, (AF2_norm),theta,y,'-.', theta, y1,'-.',theta,y2,'-.');

figure()
plot(theta, db(AF2_norm,'power'),theta,y,'-.', theta, y1,'-.',theta,y2,'-.');

AF2_norm_reduced = AF2_norm;
AF2_norm_reduced(db(AF2_norm_reduced,'p') < (rlim_low_db - 2)) = NaN; %remove nulls before polar plot
figure()
polarplot(deg2rad(theta), db(AF2_norm_reduced,'power'));
rlim([rlim_low_db 0])

%display HPBW
for i = 1:numel(AF2_norm)
    if AF2_norm(i) <= 0.5
        disp(2*theta(i))
        break
    end
end