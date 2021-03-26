close all

%%
c = physconst('LightSpeed');
Nx = 5; %# of elements in x-direction
Ny = 5; %# of elements in y-direction

I = ones(Nx,Ny); % Matri of element excitation amplitudes
Ai = ones(Nx,Ny); %Feeding Coefficient
Aiy = Aix;
Dk = 1; %substrate relative permittivity/dielectric constant
freq0 = 27e9; %Operating frequency of 28 GHz
lambda = c/freq0/sqrt(Dk);
k = 2*pi/lambda; %phase constant

% Angle domain
theta = linspace(0,pi/2,500); 
phi = linspace(0,2*pi,500);
[Theta, Phi] = meshgrid(theta, phi);

d = 1*lambda/2;   %uniform element spacing in meters
D = (0:1:Nx-1)*d;   %array of element positions
Dx = (0:1:Nx-1)'*d*ones(1,Ny); %cols=vector of distances from y-axis, rows=which column of elements
Dy = ones(Nx,1)*(0:1:Ny-1)*d; %cols=which row of elements, rows= vector of distances from x-axis

scan_theta = deg2rad(30); %scan angle of phased array
scan_phi   = deg2rad(45); %scan angle of phased array

%Cartesian axes
x_cart = sin(Theta).*cos(Phi);
y_cart = sin(Theta).*sin(Phi);

% Per elment scanning angle phase shifts
beta_x = sin(scan_theta)*cos(scan_phi);
beta_y = sin(scan_theta)*sin(scan_phi);

% Phase shifted axes used for calculation;
Gam_x = sin(Theta).*cos(Phi); 
Gam_y = sin(Theta).*sin(Phi); 

%smallest dB value shown on rad patter; used to remove extreme nulls from data 
rlim_low_db = -40; 
%% Calculate AF
AFx = zeros(size(Gam_x));
AFy = zeros(size(Gam_x));

for col= 1:size(Dx,2)
    %reshape the column of Dx into a 1-1-Nx matrix
    %each element then multiplies the 2D Gam_x, creating a 3D matrix
    %Sum along the third dimension to return a 2D matrix
    amplitude = reshape(Ai(:,col),1,1,[]);
    dx = reshape(Dx(:,col),1,1,[]);
    AFx = AFx + sum(amplitude.*exp(1i*k*dx.*(sin(Theta).*cos(Phi)-beta_x)),3);
end

for row = 1:size(Dy,1)
    amplitude = reshape(Ai(row,:),1,1,[]);
    dy = reshape(Dy(row,:),1,1,[]);
    AFy = AFy + sum(amplitude.*exp(1i*k*dy.*(sin(Theta).*sin(Phi)-beta_y)),3);
end

AF = AFx .* AFy;
AF_max = max(AF,[],'all');
AF_norm = AF./AF_max;

AF2 = abs(AF).^2;
AF2_max = max(AF2,[],'all');
AF2_norm = AF2./AF2_max;

%% Figure 1 plots spatially wihtin unit circle
figure() 
s1 = surf(sin(Theta).*cos(Phi), sin(Theta).*sin(Phi),abs(AF_norm));
s1.EdgeColor = 'none';

%% Figure 2 plots in cylindrical coordinates 
% Radius maps to elevation scanning angle
% Azimuthal angle maps to aimuthal scanning angle
% Height (z) maps to abs(AF)

x_cyl = cos(Phi).*Theta;
y_cyl = sin(Phi).*Theta;

cyl_ax_limits = [-1*pi/2 1*pi/2];
figure() % Figure 2
s2 = surf(x_cyl, y_cyl,abs(AF_norm));
s2.EdgeColor = 'none';
xlim(cyl_ax_limits)
ylim(cyl_ax_limits)
xlabel('Phi=0 Plane')
ylabel('Phi=$\pi/2$ Plane', 'Interpreter', 'latex')

%% Converting to Spherical
sph_ax_lims = [-1 1];
[x, y, z] = sph2cart(Phi,Theta,abs(AF_norm));
figure() % Figure 3
s3 = surf(x,y,z);
xlim(sph_ax_lims)
ylim(sph_ax_lims)
% zlim([0 1])
s3.EdgeColor = 'none';

figure()
patternCustom(abs(AF_norm), rad2deg(theta), rad2deg(phi), 'CoordinateSystem', 'polar')
