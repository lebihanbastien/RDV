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
nro_T = nro_interpolation(cr3bp, nro_T, nro_init_EML2, default, cst, 'Az', 70000); %Az o periselene??

%% Orbit CHASER
%Initialize NRO 
nro_C = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro_C = nro_interpolation(cr3bp, nro_C, nro_init_EML2, default, cst, 'Az', 69900);


%% Nonlinear Relative dynamics in synodic frame

%initial conditions
% theta_T = 350*pi/180;
% theta_C = 315*pi/180;

% worst case C-W
theta_T = 0*pi/180;
theta_C = 10*pi/180;

% selection of method
choice = [1 0 1 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation

% selection of hold points and TOFs
hold_points_dim = [];
TOF_dim  = [3600*6.5];
TOF_cont = [3600]; 
int_time = 300;

% other settings
toll = 1e-12;
it_max = 100;
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);


int = 1:100; %n di intervalli in cui è diviso il TOF per Luquette & CW

n = length(int);
for i=1:n
    
    intervals = int(i);
    
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

err_CW(i) = output_ch.linear.CW.err*cr3bp.L;
cc_CW(i) = output_ch.linear.CW.cc;

% err_LR(i) = output_ch.linear.LR.err*cr3bp.L;
% cc_LR(i) = output_ch.linear.LR.cc;

end


figure
plot(int,err_CW)
xlabel('No. of intervals')
ylabel('C-W dimensional error [km]')

figure
plot(int,cc_CW)
xlabel('No. of intervals')
ylabel('C-W rcond(\Phi)')

% figure
% plot(int,err_LR)
% xlabel('No. of intervals')
% ylabel('LR dimensional error [km]')
% 
% figure
% plot(int,cc_LR)
% xlabel('No. of intervals')
% ylabel('rcond(\Phi) (LR)')
