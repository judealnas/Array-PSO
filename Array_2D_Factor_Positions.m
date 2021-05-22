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
d = 1/2*lambda;   %uniform element spacing in meters
D = (0:1:Nx-1)*d;   %array of element positions
[Dx, Dy] = meshgrid(D);

% scanning angle
scan_theta_deg = 20;
scan_phi_deg = 0;
scan_theta = deg2rad(scan_theta_deg); %scan angle of phased array
scan_phi   = deg2rad(scan_phi_deg); %scan angle of phased array

%Cartesian axes
x_cart = sin(Theta).*cos(Phi);
y_cart = sin(Theta).*sin(Phi);

% Per elment scanning angle phase shifts
Beta_x = -k*Dx*sin(scan_theta)*cos(scan_phi);
Beta_y = -k*Dy*sin(scan_theta)*sin(scan_phi);
Beta = Beta_x + Beta_y;

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
        beta = Beta(j,i); % total phase of current antenna element
        AF = AF + a*exp(1i*(k*(dx*x_cart + dy*y_cart) + beta));
    end
end

AF_max = max(AF,[],'all');
AF_norm = AF./AF_max;

AF2 = abs(AF).^2;
[AF2_max] = max(AF2,[],'all');
AF2_norm = AF2./AF2_max;

%% Peak Detection
input = abs(AF_norm); % select which array factor form to process

% obtaining all peaks and their indices
TF = imregionalmax(input);  
TF_1 = islocalmax(input,1); % local maxima along dimension 1
TF_2 = islocalmax(input,2); % local maxima along dimension 2
TF_intersect = TF_1 & TF_2; % Maxima that exist in both dimensions

% vector subset of just the detected local maxima amplitudes and their
% linear indices relative to input
mag_local_maxima = input(TF_intersect); 
idx_local_maxima = find(TF_intersect);

% searches subset of local maxima for the global maxima and its linear
% index RELATIVE TO SUBSET. Use this to index idx_local_maxima to obtain
% the linear index of the global max relative to input
[glob_max, tempidx_glob_max] = max(mag_local_maxima);
idx_glob_max = idx_local_maxima(tempidx_glob_max);

% Remove global max index from idx_local_maxima. 
% Largest of remaining peaks is SLL.
idx_side_maxima = idx_local_maxima(idx_local_maxima ~= idx_glob_max);
[SLL, tempidx_SLL] = max(input(idx_side_maxima));
idx_SLL = idx_local_maxima(tempidx_SLL);

%% Plots against angle in rectangular
figure();
patternCustom(input, rad2deg(theta), rad2deg(phi), 'CoordinateSystem', 'rectangular')
% s0 = surf(Theta, Phi, input);
% s0.EdgeColor = 'none';
xlabel('Theta', 'Interpreter', 'latex')
ylabel('Phi', 'Interpreter', 'latex')
hold on

% plot points at all detected maxima
plot3(rad2deg(Theta(TF_1)),...
      rad2deg(Phi(TF_1)),...
      input(TF_1), 'y.', 'MarkerSize', 15);

% plot points at all detected maxima
plot3(rad2deg(Theta(TF_2)),...
      rad2deg(Phi(TF_2)),...
      input(TF_2), 'y.', 'MarkerSize', 15);

% plot points at all detected maxima
plot3(rad2deg(Theta(TF_intersect)),...
      rad2deg(Phi(TF_intersect)),...
      input(TF_intersect), 'y.', 'MarkerSize', 15);

% points at the global max
plot3(rad2deg(Theta(idx_glob_max)),...
      rad2deg(Phi(idx_glob_max)),...
      input(idx_glob_max), 'g.','MarkerSize', 20);
  
% points at detected SLL
plot3(rad2deg(Theta(idx_SLL)),...
      rad2deg(Phi(idx_SLL)),...
      input(idx_SLL), 'rd','MarkerSize', 25);
hold off
%% Figure plots spatially wihtin unit circle
figure() 
z_cart = input;
s1 = surf(x_cart, y_cart,z_cart);
s1.EdgeColor = 'none';
xlabel("X");
ylabel("Y");
hold on

% local maxima in dimension 1
plot3(x_cart(TF_1),...
      y_cart(TF_1),...
      z_cart(TF_1), 'k.', 'MarkerSize', 15);

% local maxima in dimension 2
plot3(x_cart(TF_2),...
      y_cart(TF_2),...
      z_cart(TF_2), 'b.', 'MarkerSize', 15);

% points at intersection
plot3(x_cart(TF_intersect),...
      y_cart(TF_intersect),...
      z_cart(TF_intersect), 'r.', 'MarkerSize', 50);

% points at the global max
plot3(x_cart(idx_glob_max),...
      y_cart(idx_glob_max),...
      z_cart(idx_glob_max), 'g.','MarkerSize', 20);
  
% points at detected SLL
plot3(x_cart(idx_SLL),...
      y_cart(idx_SLL),...
      z_cart(idx_SLL), 'rd','MarkerSize', 15);

% plot3(x_cart(TF_intersect), y_cart(TF_intersect), z_cart(TF_intersect), 'r.','MarkerSize', 25);


hold off
%% Plots in cylindrical coordinates 
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

%% Plotting in Spherical
figure()
patternCustom(input, rad2deg(theta), rad2deg(phi), 'CoordinateSystem', 'polar')
