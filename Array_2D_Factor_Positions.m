close all

%%
c = physconst('LightSpeed');
Nx = 5; %# of elements in x-direction
Ny = 5; %# of elements in y-direction

I = ones(Ny,Nx); % Matri of element excitation amplitudes
Ai = ones(Ny,Nx); %Feeding Coefficient

% Constants associated with material
Dk = 1; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant

% Angle domain; covers top hemisphere
theta = linspace(0,pi/2,500); 
phi = linspace(0,2*pi,500);
[Theta, Phi] = meshgrid(theta, phi);

% spacing matrices
% Dx: size = [Nx Ny]; The value in the ith row, jth column
%            is the distance in x from the y-axis of the element
%            in the corresponding row & column
% Dy: size = [Nx Ny]; The value in the ith row, jth column
%            is the distance in y from the x-axis of the element
%            in the corresponding row & column 
d = 1*lambda;   %uniform element spacing in meters
D = (0:1:Nx-1)*d;   %array of element positions
[Dx, Dy] = meshgrid(D);

% scanning angle
scan_theta = deg2rad(0); %scan angle of phased array
scan_phi   = deg2rad(-50); %scan angle of phased array

%Cartesian axes
x_cart = sin(Theta).*cos(Phi);
y_cart = sin(Theta).*sin(Phi);

% Per elment scanning angle phase shifts
beta_x = sin(scan_theta)*cos(scan_phi);
beta_y = sin(scan_theta)*sin(scan_phi);
Beta_x = -k*Dx*sin(scan_theta)*cos(scan_phi);
Beta_y = -k*Dy*sin(scan_theta)*sin(scan_phi);

%smallest dB value shown on rad patter; used to remove extreme nulls from data 
rlim_low_db = -40; 
%% Calculate AF
AF = zeros(size(x_cart));

% iterate over each of the Nx*Ny elements, calculating
% and adding the 2D array factor contributions
for i = 1:Nx % iterate in x-direction (i.e. cols or dim 2)
    for j = 1:Ny %iterate in y-direction (i.e. rows or dim 1)
        dx = Dx(j,i); % x location of current antenna element
        dy = Dy(j,i); % y location of current antenna element
        a = Ai(j,i);  % feed amplitude of current antenna element
        beta = Beta_x(j,i) + Beta_y(j,i); % total phase of current antenna element
        AF = AF + a*exp(1i*(k*(dx*x_cart + dy*y_cart) + beta));
    end
end

AF_max = max(AF,[],'all');
AF_norm = AF./AF_max;

AF2 = abs(AF).^2;
[AF2_max] = max(AF2,[],'all');
AF2_norm = AF2./AF2_max;

%% Peak Detection
input = abs(AF2_norm); % select which array factor form to process

% obtaining all peaks and their indices
TF = imregionalmax(input);
[ind_max_row, ind_max_col] = find(TF);
ind_peaks = [ind1, ind2];

% removing the global peak
peaks = input(TF);
subpeaks = find(peaks(peaks < max(peaks));
[SLL, ind_SLL]= max(peaks(peaks ~= input_max));
ind_SLL = [ind_max_input(ind_SLL,1) ind_max_input(ind_SLL,2)];
%% Figure 1 plots spatially wihtin unit circle
figure() 
z_cart = input;
s1 = surf(x_cart, y_cart,z_cart);
s1.EdgeColor = 'none';
hold
plot3(x_cart(TF), y_cart(TF), z_cart(TF), 'r+'); % plotting peak locations determined by imregionmax(0
plot3(x_cart(ind_max_input(1),ind_max_input(2)),...
      y_cart(ind_max_input(1),ind_max_input(2)),... 
      z_cart(ind_max_input(1),ind_max_input(2)), 'gs'); % plotting the max location
plot3(x_cart(ind_SLL(1),ind_SLL(2)),...
      y_cart(ind_SLL(1),ind_SLL(2)),...
      z_cart(ind_SLL(1),ind_SLL(2)), 'bp'); % plotting second highest peaks
  hold
%% Figure 2 plots in cylindrical coordinates 
% Radius maps to elevation scanning angle
% Azimuthal angle maps to aimuthal scanning angle
% Height (z) maps to abs(AF)

x_cyl = cos(Phi).*Theta;
y_cyl = sin(Phi).*Theta;
z_cyl = input;
cyl_ax_limits = [-1*pi/2 1*pi/2];
figure() % Figure 2
s2 = surf(x_cyl, y_cyl,z_cyl);
s2.EdgeColor = 'none';
xlim(cyl_ax_limits)
ylim(cyl_ax_limits)
xlabel('Phi=0 Plane', 'Interpreter', 'latex')
ylabel('Phi=$\pi/2$ Plane', 'Interpreter', 'latex')

%% Converting to Spherical
% sph_ax_lims = [-1 1];
% [x, y, z] = sph2cart(Phi,Theta,input);
% figure() % Figure 3
% s3 = surf(x,y,z);
% xlim(sph_ax_lims)
% ylim(sph_ax_lims)
% % zlim([0 1])
% s3.EdgeColor = 'none';

figure()
patternCustom(input, rad2deg(theta), rad2deg(phi), 'CoordinateSystem', 'polar')
