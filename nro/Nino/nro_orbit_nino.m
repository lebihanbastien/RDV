%--------------------------------------------------------------------------
% Plots an NRO orbit.
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = true;  %plot also the results in X-Y plane
default.plot.XZ          = false;  %plot also the results in X-Z plane
default.plot.YZ          = false;  %plot also the results in Y-Z plane
default.plot.TD          = true;  %plot also the results in three-D


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

% initial conditions
% theta_T = 350*pi/180;
% theta_C = 315*pi/180;

theta_T = 180*pi/180;
theta_C = 178*pi/180;

% selection of method
choice = [1 0 0 0]; %[1 0 0 0] CW, [0 1 0 0] SL, [0 0 1 0] LR, [1 1 1 0] previous 3, [0 0 0 1] continuation

% selection of hold points and TOFs
hold_points_dim = [-5 0 -1];
TOF_dim  = [3600*2 3600];
TOF_cont = [3600 600]; 
int_time = 300;

% other settings
intervals = 100; %n di intervalli in cui � diviso il TOF per Luquette & CW
toll = 1e-12;
it_max = 100;
options = odeset('Reltol', default.ode113.RelTol, 'Abstol', default.ode113.AbsTol);

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

output_ch = Rendezvous_choice(choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, 1);

%output_ch = Rendezvous_choice(hold_points_dim, TOF_dim, theta_T, theta_C, nro_T, nro_C, cr3bp, toll, it_max, intervals, TOF_cont, int_time, options, choice, 1);
%output = Rendezvous(hold_points_dim, TOF_dim, t_T, t_C, nro_T, nro_C, cr3bp, toll, it_max, intervals, options, 1);
%output = Rendezvous_LVLH (hold_points_dim, TOF_dim, t_T, t_C, nro_T, nro_C, cr3bp, toll, it_max, options, 1);

%% Plotting Results

% %Lambert solution in LVLH
% np = size(hold_points_dim,1) + 2;
% figure
%  for i=1:np-1
%      p1 = plot3(output.Lambert_LVLH(i).yv(:,1),output.Lambert_LVLH(i).yv(:,2),output.Lambert_LVLH(i).yv(:,3),'r');
%      hold on 
%      axis equal
%      grid on
%      p2 = plot3(output2.Lambert_LVLH(i).transfer.yv(:,1),output2.Lambert_LVLH(i).transfer.yv(:,2),output2.Lambert_LVLH(i).transfer.yv(:,3),'go');
%      p4(i) = plot3(output.Lambert_LVLH(i).yv(end,1),output.Lambert_LVLH(i).yv(end,2),output.Lambert_LVLH(i).yv(end,3),'bo');
%  end
%   
%  p3 = plot3(output.Lambert_LVLH(1).yv(1,1),output.Lambert_LVLH(1).yv(1,2),output.Lambert_LVLH(1).yv(1,3),'co');
%  p5 = plot3(0,0,0,'ko');
%  
%  if np>2  
%    for i=1:np-2
%       str(i) = cellstr(sprintf('%s %d','HP',i));
%    end 
%    pp = [p1 p2 p3 p4(1:np-2) p5];
%    nom = ['Lambert with Syn Dynamics' 'Lambert with LVLH Dynamics' 'Start' str 'Docking'];
%    legend(pp,nom)
% else
%    legend([p1 p2 p3 p5],'Lambert with Syn Dynamics', 'Lambert with LVLH Dynamics', 'Start','Docking')
% end
%  
%  title('Lambert Arcs in LVLH Frame')
%  xlabel('x_{LVLH}')
%  ylabel('y_{LVLH}')
%  zlabel('z_{LVLH}')
 


%% Dimensional Error for both dynamics

%  figure
%  for i=1:np-1
%  semilogy(i,output.CW.Lambert_syn.error(i).tr(1).*cr3bp.L,'ro')
%  hold on
%  grid on
%  semilogy(i,output2.CW.Lambert_LVLH.error(i).tr(1).*cr3bp.L,'rd')
%  semilogy(i,output.SL.Lambert_syn.error(i).tr(1).*cr3bp.L,'bo')
%  semilogy(i,output2.SL.Lambert_LVLH.error(i).tr(1).*cr3bp.L,'bd')
%  semilogy(i,output.LR.Lambert_syn.error(i).tr(1).*cr3bp.L,'go')
%  semilogy(i,output2.LR.Lambert_LVLH.error(i).tr(1).*cr3bp.L,'gd')
%  end
%  
%  set(gca,'XTick',[1:i] );
%  legend('CW_{syn}','CW_{LVLH}','SL_{syn}','SL_{LVLH}','LR_{syn}','LR_{LVLH}')
%  title('Initial Dimensional Error for each Transfer')
%  
%  xlabel('Transfer')
%  ylabel('Error [km]')

 
