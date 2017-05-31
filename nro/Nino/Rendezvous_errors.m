function output = Rendezvous_errors (choice, hold_points_dim, TOF_dim, init_cond, orbits, settings, continuation, plot_yes)
% Author: Antonino Campolo 2016/2017 (with the Holy help of the Fabolous Mr. BLB) 

% This code computes the rendezvous trajectory for a user-given number of
% hold points, initial orbits and TOFs 

% hold_points_dim is a matrix whose rows are the LVLH DIMENSIONAL (km) positions of the hold
% points. If no hold points are needed, hold_points_dim =[];

% TOF_dim is a vector which contains the (n. of hold points +1) DIMENSIONAL (s) time of flights for
% each transfer from one hold point to the next. At least one TOF is needed
% to compute the rendezvous

% theta_T & theta_C (radians) are the "true anomalies" of the target and the chaser
% on the planar orbit resulting from the projection of the 3D chosen orbit
% on the plane formed by its position and velocity vector at the periselene

% nro_T & nro_C are the structures regarding target and chaser's orbits
% (from main)

% toll is the tolerance on the error in the Lambert cycle. The error is
% computed as the norm of the vectorial difference between the integrated
% final position and the desired final position

% it_max is the maximum number of iteration allowed for the Lambert's cycle

% intervals is the number of intervals in which the Luquette's transfer is
% divided. For each interval, the STM is evaluated in order to compose a
% final STM which is better conditioned

% TOF_cont is a vector which contains the starting DIMENSIONAL (s) Time of 
% Flights for each transfer for the continuation method. 

% int_time is the DIMENSIONAL (s) time separation between two evaluation of
% Lambert's arcs in the continuation algorithm

% choice is a 3x1 vector wich selects a particular first guess for the
% Lambert's cycle. If 
% choice = [1 0 0 0], CW is selected.
% choice = [0 1 0 0], SL is selected
% choice = [0 0 1 0], LR is selected
% choice = [1 1 1 0], all three previous methos are used to generate an
% hybrid solution
% choice = [0 0 0 1], continuation algorithm is selected
% No other combinations are possible

% if plot is required, plot_yes = 1. Otherwise, put any other number



%% Inputs generation

theta_T = init_cond.T;
theta_C = init_cond.C;
TOF_cont = continuation.TOF_cont;
int_time = continuation.int;
nro_T = orbits.T;
nro_C = orbits.C;
it_max = settings.it_max;
cr3bp = settings.cr3bp;


% if isempty(choice) || length(nonzeros(choice)) == 2
%    error('Select one or all the first guesses')
% end

t_T = interp1(nro_T.alpha,nro_T.tv,theta_T,'spline');
t_C = interp1(nro_C.alpha,nro_C.tv,theta_C,'spline');

% BLB dovrebbe mostrare che in realtà il tutto non dipende da t.
t      = t_T;                     %defines the initial relative position of synodic/inertial frame
x0_T   = state_time(t_T,nro_T)';
x0_C   = state_time(t_C,nro_C)';
x0_rel = x0_C - x0_T;             %synodic frame 
phi0   = reshape(eye(6,6),[],1);
 
% states in LVLH
x0_LVLH     = rsyn2rlvlh(t, x0_rel, x0_T, cr3bp.mu);   %adim
hold_points = hold_points_dim./cr3bp.L;                %LVLH
dock        = [0 0 0];

%fof = x0_LVLH(1:3)*cr3bp.L
% successive positions and velocities in the LVLH frame
r_f = [x0_LVLH(1:3)'; hold_points; dock];    
np  = size(r_f,1);
v_f = [x0_LVLH(4:6)'; zeros(np-1,3)];  
 
% TOF and times vector
TOF = TOF_dim.*2*pi/cr3bp.T;


% possible errors
if length(TOF) ~= np-1
   error('The number of TOFs must be equal to the n.of hold points + 1')
end

% if length(TOF_cont) ~= np-1
%    error('The number of TOFs for the continuation method must be equal to the n.of hold points + 1')
% end

% for i=1:np-1
%     if TOF_cont(i).*2*pi/cr3bp.T >= TOF(i)
%        error('TOF_cont must be < than TOF for each transfer')
%     end
% end

% times for each hold point
t_f    = zeros(np,1);
t_f(1) = t;   %penso si possa mettere direttamente t_T. Aspetta che BLB dimostri la non dip da t.
for i=1:np-1
    t_f(i+1) = t+sum(TOF(1:i));  %t must be adimensional
end

% states in Synodic frame
m_hp   = size(hold_points,1);
hp_syn = zeros(m_hp,6);

if m_hp > 0
   for i=1:m_hp
       hp_syn(i,:) = rlvlh2rsyn(t_f(i+1), [hold_points(i,:) 0 0 0]', state_time(t_f(i+1),nro_T)', cr3bp.mu);
   end
   
   r_f_syn = [x0_rel(1:3)'; hp_syn(:,1:3); dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)]; 
 
else
   hp_syn = []; 
   r_f_syn = [x0_rel(1:3)'; hp_syn; dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-1,3)]; 
end
 

%% First guesses & Lambert
for i = 1:np-1
    r_LVLH = r_f(i:i+1,:);
    v_LVLH = v_f(i,:);
    r_syn  = r_f_syn(i:i+1,:);
    r_next = r_syn(2,:);
    v_syn  = v_f_syn(i,:);
    TOF_ca = TOF_cont(i);
    T      = TOF(i);
    tt     = t_f(i);
    tt_2   = t_f(i+1);
    
    LVLH.r         = r_LVLH;
    LVLH.v         = v_LVLH;
    SYN.r          = r_syn; 
    SYN.v          = v_syn;
    SYN.r_next     = r_next;
    times_lin.T    = T; 
    times_lin.tt   = tt;
    times_lin.tt_2 = tt_2;
    times_cont.TOF = TOF_ca;
    times_cont.int = int_time;
    
    pop = Terminator(choice, LVLH, SYN, times_lin, times_cont, settings, nro_T, phi0);
    sol(i) = pop;
end

 %% Building Hybrid Solution
%Lam = building(choice, sol, np, it_max);
 
%% deltaV calculation
%dV = dV_eval(Lam, x0_rel, np, cr3bp);

%% Outputs
output.linear = sol;
output.init_LVLH = x0_LVLH;
%output.Lambert = Lam;

% deltaV
%output.dV = dV;
 
%% Plotting

if plot_yes == 1
   %plot_1   %linear models
    plot_3_choice(choice, sol, Lam, r_f_syn, np, cr3bp, nro_T, nro_C, x0_T, x0_C)  %Lambert
end

 %la domanda da porsi è: utilizzo le eqz adim o dim? Privilegio la rapidità
 %del codice o la precisione numerica? Rispondere dopo. Adesso
 %privilegiare l'omogeneità del codice.
 