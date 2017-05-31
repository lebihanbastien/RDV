%--------------------------------------------------------------------------
% Plots an NRO orbit.
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = false;  %plot also the results in X-Y plane
default.plot.XZ          = false;  %plot also the results in X-Z plane
default.plot.YZ          = false;  %plot also the results in Y-Z plane
default.plot.TD          = false;   %plot also the results in three-D
%% Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit TARGET
%Initialize NRO 
nro_T = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro_T = nro_interpolation(cr3bp, nro_T, nro_init_EML2, default, cst, 'altitudeOfPerigee', 500); %Az o periselene??

%% Orbit CHASER
%Initialize NRO 
nro_C = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
% nro_C = nro_interpolation(cr3bp, nro_C, nro_init_EML2, default, cst, 'altitudeOfPerigee', 499.5);
 nro_C = nro_interpolation(cr3bp, nro_C, nro_init_EML2, default, cst, 'altitudeOfPerigee', 499.74);


%% Nonlinear Relative dynamics in synodic frame

%initial conditions
% theta_T = 350*pi/180;
% theta_C = 315*pi/180;

% selection of method
choice = [1 1 1 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation

% selection of hold points and TOFs
hold_points_dim = [];
TOF_cont = 3600*2; 
int_time = 300;

% other settings
intervals = 100; %n di intervalli in cui è diviso il TOF per Luquette
toll = 1e-12;
it_max = 100;
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

% inputs for cycles (change orbits for different points, rel dist from 500
% m)
targ = 90; %target angular position in degrees
step = 0.01; %step for angular cycle
dtheta = 3.57; % +- targ position
p_m = targ - dtheta;
p_p = targ + dtheta;
if p_m*p_p < 0 
   p_m = 360 - abs(p_m);
   CC = [p_m:step:360-step,0:step:p_p]; %non vero in generale
else
CC = p_m:step:p_p;


end
TT = 600:300:36000;
% TT = 600:35400:36000;






% target position (fixed)
theta_T = targ*pi/180;

% initializing
m = length(CC);
n = length(TT);

err_CW = zeros(m,n);
err_SL = zeros(m,n);
err_LR = zeros(m,n);
x_init = zeros(m,n);
y_init = zeros(m,n);
z_init = zeros(m,n);
init_dist = zeros(m,n);


for i = 1:m
    for j = 1:n
        theta_C  = CC(i)*pi/180;
        TOF_dim  = TT(j);
        
        CC(i)
        TT(j)

%structures as inputs
init_cond.T = theta_T;
init_cond.C = theta_C;
continuation.TOF_cont = TOF_cont;
continuation.int = int_time;
orbits.T = nro_T;
orbits.C = nro_C;
settings.Luq_int = intervals;
settings.toll = toll;
settings.it_max = it_max;
settings.options = options;
settings.cr3bp = cr3bp;

output_ch = Rendezvous_errors(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, 0);

err_CW(i,j)    = output_ch.linear.CW.err*cr3bp.L;
err_SL(i,j)    = output_ch.linear.SL.err*cr3bp.L;
err_LR(i,j)    = output_ch.linear.LR.err*cr3bp.L;
x_init(i,j)    = output_ch.init_LVLH(1)*cr3bp.L;
y_init(i,j)    = output_ch.init_LVLH(2)*cr3bp.L;
z_init(i,j)    = output_ch.init_LVLH(3)*cr3bp.L;
init_pos       = [x_init(i,j) y_init(i,j) z_init(i,j)];
init_dist(i,j) = norm(init_pos);

    end
end

%% plot
CC_plot = -dtheta:step:dtheta;
TT_plot = TT/3600;


figure
surf(CC_plot,TT_plot,err_CW','edgecolor', 'none');
colorbar
%view(0,0)
xlabel('\Delta \theta [°]')
ylabel('TOF [h]')
zlabel('C-W dimensional error [km]')

% hold on
% % [r,c] = find(err_CW == max(max(err_CW)))
% r = 201;
% c = 77;
% plot3(CC_plot(r), TT_plot(c), err_CW(r,c)','ro','LineWidth',2)

figure
surf(CC_plot,TT_plot,err_SL','edgecolor', 'none');
colorbar
%view(0,0)
xlabel('\Delta \theta [°]')
ylabel('TOF [h]')
zlabel('SL dimensional error [km]')

figure
surf(CC_plot,TT_plot,err_LR','edgecolor', 'none');
colorbar
%view(0,0)
xlabel('\Delta \theta [°]')
ylabel('TOF [h]')
zlabel('LR dimensional error [km]')

% figure
% surf(CC_plot,TT_plot,init_dist','edgecolor', 'none')
% colorbar
% xlabel('\Delta \theta [°]')
% ylabel('TOF [h]')
% zlabel('Initial relative distance[km]')
% 
% a = min(min(init_dist))
% b = max(max(init_dist))
