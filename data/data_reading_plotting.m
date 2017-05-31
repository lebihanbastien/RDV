%--------------------------------------------------------------------------
% This matlab file reads the data in txt format and outputs them in matlab
% format. Additionnal plots are produced to check the monotony of the Az =
% f(x0) function.
% RM: the data files contain the NORTHERN family.
%
% Author: BLB
% Version: 1.0
% Year: 2014
%--------------------------------------------------------------------------
% Add subfolder to the path
addpath('./data');
addpath('./computation');
addpath('./init');

%% Constants init
cst = constants_init();

%% Parameters init (to default values, see values within routine)
default = parameters_default_init(cst);

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% FILE reading (L2)
fileID = fopen('halo_init_matrix_EML2.txt','r');
C_data1 = textscan(fileID,'%f %f %f %f %f %f %f','CollectOutput',1);
fclose(fileID);
halo_init_EML2.matrix = C_data1{1};
halo_init_EML2.Azlimit = 0.19953;  %ensure monotony of Az = f(x0) and good interpolation

%Append jacobi vector
for i=1:size(halo_init_EML2.matrix,1)
   halo_init_EML2.matrix(i,8) = jacobi(halo_init_EML2.matrix(i,1:6), cr3bp.mu);
end
% Limits to ensure good interpolation
halo_init_EML2.Cjaclimit(1) = 3.0320;
halo_init_EML2.Cjaclimit(2) = 3.1641;

% Save
save halo_init_matrix_EML2 halo_init_EML2

%% FILE reading (L1)
fileID = fopen('halo_init_matrix_EML1.txt','r');
C_data1 = textscan(fileID,'%f %f %f %f %f %f %f','CollectOutput',1);
fclose(fileID);
halo_init_EML1.matrix = C_data1{1};
halo_init_EML1.Azlimit = 0.32983;  %ensure monotony of Az = f(x0) and good interpolation

%Append jacobi vector
for i=1:size(halo_init_EML1.matrix,1)
   halo_init_EML1.matrix(i,8) = jacobi(halo_init_EML1.matrix(i,1:6), cr3bp.mu);
end
% Limits to ensure good interpolation
halo_init_EML1.Cjaclimit(1) = 2.9213;
halo_init_EML1.Cjaclimit(2) = 3.1862;


save halo_init_matrix_EML1 halo_init_EML1

%% Plot (L2)
figure;
hold on;
grid on;
title('Initial condition for halo orbit computation @EML2: Az as a function of x0', 'FontSize', 20);
xlabel('x0 [km]', 'FontSize', 20);
ylabel('Az [km]', 'FontSize', 20);
plot(halo_init_EML2.matrix(:,1), halo_init_EML2.matrix(:,7), 'LineWidth', 3);
set(gca,  'FontSize', 20);


figure;
hold on;
grid on;
title('Jacobi constant as a function of x0 @EML2', 'FontSize', 20);
ylabel('Jacobi [-]', 'FontSize', 20);
xlabel('x0 [km]', 'FontSize', 20);
plot(halo_init_EML2.matrix(:,1), halo_init_EML2.matrix(:,8), 'LineWidth', 3);
set(gca,  'FontSize', 20);
%% Plot (L1)
figure;
hold on;
grid on;
title('Initial condition for halo orbit computation @EML1: Az as a function of x0', 'FontSize', 20);
xlabel('x0 [km]', 'FontSize', 20);
ylabel('Az [km]', 'FontSize', 20);
plot(halo_init_EML1.matrix(:,1), halo_init_EML1.matrix(:,7), 'LineWidth', 3);
set(gca,  'FontSize', 20);

figure;
hold on;
grid on;
title('Jacobi constant as a function of x0 @EML1', 'FontSize', 20);
ylabel('Jacobi [-]', 'FontSize', 20);
xlabel('x0 [km]', 'FontSize', 20);
plot(halo_init_EML1.matrix(:,1), halo_init_EML1.matrix(:,8), 'LineWidth', 3);
set(gca,  'FontSize', 20);