close all

%%
c = physconst('LightSpeed');
Nx = 5; %# of elements in x-direction
Ny = 5; %# of elements in y-direction

Ai = ones(Nx,Ny); %Feeding Coefficient
Dk = 1; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant

theta = (-90:1:90); %degree angle
phi = (0:1:360);

Gam_x = cos(deg2rad(phi(:)))*sin(deg2rad(theta(:)))';
Gam_y = sin(deg2rad(phi(:)))*sin(deg2rad(theta(:)))';

d = 2*lambda/2;   %uniform element spacing in meters

D = (0:1:Nx-1)*d;   %array of element positions
Dx = (0:1:Nx-1)'*d*ones(1,Ny); %cols=vector of distances from y-axis, rows=which column of elements
Dy = ones(Nx,1)*(0:1:Ny-1)*d; %cols=which row of elements, rows= vector of distances from x-axis

scan_theta_deg = 30; %scan angle of phased array
scan_phi_deg = 30; %scan angle of phased array
rlim_low_db = -40; %smallest dB value shown on rad patter; used to remove extreme nulls from data 
%% Calculate AF
AFx = zeros(size(Gam_x));
AFy = zeros(size(Gam_x));

for col= 1:size(Dx,2)
    %reshape the column of Dx into a 1-1-Nx matrix
    %each element then multiplies the 2D Gam_x, creating a 3D matrix
    %Sum along the third dimension to return a 2D matrix
    AFx = AFx + sum(exp(1i*k*reshape(Dx(:,col),1,1,[]).*Gam_x),3);
end

for row = 1:size(Dy,1)
    AFy = AFy + sum(exp(1i*k*reshape(Dy(row,:),1,1,[]).*Gam_y),3);
end

AF = Ai .* AFx .* AFy;
AF2 = abs(AF/max(AF)).^2;
AF2_norm = AF2./max(AF2);

s = surf(Gam_x, Gam_y, AF2);
s.EdgeColor = 'none';
%% Peak Finding
%Locate indices of peaks; remove index of main lobe; calculate max of remaining peaks

%max SLL using valley finder; valleys used to find main node null points
% ind_vall = find(islocalmin(AF2)); %find indices of valleys 
% n = sort([ind_vall(:); ind_max(:)]); %sorted array of indices of valleys and the main peaks
% null1 = n(find(n == ind_max(1)) + 1); %ind_max(1) = index of boresight
% null2 = n(find(n == ind_max(2)) - 1); %ind_max(2) = index of rear boresight
% 
% [max_SLL1,ind_max_SLL] = max(AF2(null1:null2))
% 
% %max SLL using peak finder
% ind_peak = find(islocalmax(AF2)); %find indices of peaks (detection may cause problems)
% 
% ind_side = setdiff(ind_peak,ind_max); %remove index of boresights; gives indices of side lobe peaks
% 
% [max_SLL2, ind2_max_SLL] = max(AF2(ind_side)); %ind2 = second-order index, an index to retrieve ANOTHER index
% max_SLL2 == AF2(ind_side(ind2_max_SLL));
% 
% %% Plot
% %data to plot located SL peaks
% y = zeros(size(theta));
% y(ind_side) = 1;
% 
% y1 = zeros(size(theta));
% y1(ind_peak) = 1;
% 
% y2 = zeros(size(theta));
% y2(ind_vall) = 1;
% 
% %Replacing nulls with NaN to eliminate spikes
% AF2_norm = AF2/AF2_max;
% 
% figure()
% plot(theta, (AF2_norm),theta,y,'-.', theta, y1,'-.',theta,y2,'-.');
% 
% figure()
% plot(theta, db(AF2_norm,'power'),theta,y,'-.', theta, y1,'-.',theta,y2,'-.');
% 
% AF2_norm_reduced = AF2_norm;
% AF2_norm_reduced(db(AF2_norm_reduced,'p') < (rlim_low_db - 2)) = NaN; %remove nulls before polar plot
% figure()
% polarplot(deg2rad(theta), db(AF2_norm_reduced,'power'));
% rlim([rlim_low_db 0])
% 
% %display HPBW
% for i = 1:numel(AF2_norm)
%     if AF2_norm(i) <= 0.5
%         disp(2*theta(i))
%         break
%     end
% end