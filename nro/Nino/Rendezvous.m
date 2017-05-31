function output = Rendezvous (hold_points_dim, TOF_dim, t_T, t_C, nro_T, nro_C, cr3bp, toll, it_max, intervals, options, plot_yes)

% Author: Antonino Campolo 2016 (with the Holy help of the fabolous Mr. BLB) 

% This code computes the rendezvous trajectory for a user-given number of
% hold points, initial orbits and TOFs 

% hold_points_dim is a matrix whose rows are the LVLH DIMENSIONAL (km) positions of the hold
% points. If no hold points are needed, hold_points_dim =[];

% TOF_dim is a vector which contains the (n. of hold points +1) DIMENSIONAL (s) time of flights for
% each transfer from one hold point to the next. At least one TOF is needed
% to compute the rendezvous

% t_T & t_C are ADIMENSIONAL times that PARAMETRISE the initial state vectors for target and
% chaser. t_T = t_C = 0, both T and C are at the periselene

% nro_T & nro_C are the structures regarding target and chaser's orbits
% (from main)

% toll is the tolerance on the error in the Lambert cycle. The error is
% computed as the norm of the vectorial difference between the integrated
% final position and the desired final position

% if plot is required, plot_yes==1.  the dimensiOtherwise, put any other number

%% Inputs generation

%BLB dovrebbe mostrare che in realtà il tutto non dipende da t.
t = t_T;            %defines the initial relative position of synodic/inertial frame
x0_T = state_time(t_T,nro_T)';
x0_C = state_time(t_C,nro_C)';
x0_rel = x0_C - x0_T;       %synodic frame 
phi0 = reshape(eye(6,6),[],1);
 
%states in LVLH
x0_LVLH = rsyn2rlvlh(t, x0_rel, x0_T, cr3bp.mu);   %adim
hold_points = hold_points_dim./cr3bp.L;    %LVLH
dock = [0 0 0];
 
r_f = [x0_LVLH(1:3)'; hold_points; dock]  %Docking LVLH
np = size(r_f,1);
v_f = [x0_LVLH(4:6)'; zeros(np-2,3)];   %dal docking non parte nulla (np-2)
 
%TOF and times vector
TOF = TOF_dim.*2*pi/cr3bp.T;
 
if length(TOF) ~= np-1
   error('The number of TOFs must be equal to the n.of hold points + 1')
end
 
t_f(1) = t;   %penso si possa mettere direttamente t_T. Aspetta che BLB dimostri la non dip da t.
for i=1:np-1
    t_f(i+1) = t+sum(TOF(1:i));  %t must be adimensional
end

%states in SYN
if size(hold_points,1)>0
 for i=1:size(hold_points,1)
     hp_syn(i,:) = rlvlh2rsyn(t_f(i+1), [hold_points(i,:) 0 0 0]', state_time(t_f(i+1),nro_T)', cr3bp.mu);
 end
    r_f_syn = [x0_rel(1:3)'; hp_syn(:,1:3); dock];
    v_f_syn = [x0_rel(4:6)'; zeros(np-2,3)]; 
 
else
   hp_syn = []; 
   r_f_syn = [x0_rel(1:3)'; hp_syn; dock];
   v_f_syn = [x0_rel(4:6)'; zeros(np-2,3)]; 
end
 
%attenzione al fatto che LVLH e syn_rel sono centrate sul target, quindi
%è vero che [0 0 0]_rsyn = [0 0 0]_LVLH. Tuttavia questo non è vero in generale

%% Clohessy Wiltshire: First guesses & Lambert

tic()
%Clohessy Wiltshire (LVLH)
for i=1:np-1
   [x_CW_LVLH(i,:), deltaV_CW(i), state_CW(i).yv] = Docking(r_f(i:i+1,:), v_f(i,:), TOF(i), state_time(t_f(i),nro_T), cr3bp);
    x_CW_syn(i,:) = rlvlh2rsyn(t_f(i), x_CW_LVLH(i,:)', state_time(t_f(i),nro_T)', cr3bp.mu)';
end    
x_CW_syn(size(x_CW_syn,1) + 1, :) = zeros(1,6); %I need the final position
 
%Lambert CW (SYN)
for i=1:np-1
   [Transfer_CW(i).it, it_CW(i), err_CW(i).tr] = Lambert_cycle(x_CW_syn(i:i+1,:), state_time(t_f(i),nro_T), phi0, TOF(i), toll, it_max, cr3bp, options);
end
toc()

%% Straight Line Targeting (LVLH)

tic()
for i=1:np-1
   [x_SL_LVLH(i,:), deltaV_SL(i), state_SL(i).yv] = Docking_SL(r_f(i:i+1,:), v_f(i,:), TOF(i));
    x_SL_syn(i,:) = rlvlh2rsyn(t_f(i), x_SL_LVLH(i,:)', state_time(t_f(i),nro_T)', cr3bp.mu)';
end
x_SL_syn(size(x_SL_syn,1) + 1, :) = zeros(1,6); 

 %Lambert SL (SYN)
for i=1:np-1
   [Transfer_SL(i).it, it_SL(i), err_SL(i).tr] = Lambert_cycle(x_SL_syn(i:i+1,:), state_time(t_f(i),nro_T), phi0, TOF(i), toll, it_max, cr3bp, options);
end
toc()

%% Linerarized Relative Targeting (SYN)

% for i=1:np-1
%    [x_LR_syn(i,:), deltaV_LR(i), state_LR(i).yv] = Docking_LRT(r_f_syn(i:i+1,:), v_f_syn(i,:), TOF(i), state_time(t_f(i),nro_T), cr3bp);
% end

%state_LR è complicato e pesante come output. dovrei calcolare A per ogni
%intervallino (non esiste A globale come prodotto delle singole A) e poi
%usare dei punti intermedi per ogni transfer. Diventa come un multiple
%shooting e non avrei il dV iniziale per il Lambert xk il dV iniziale
%sarebbe quello calcolato per arrivare al prossimo patch point. Visto che
%state_LR serviva solo per il plot, se serve c'è la vecchia versione
%Docking_LRT

tic()
for i=1:np-1
   [x_LR_syn(i,:), deltaV_LR(i)] = Docking_Luquette(r_f_syn(i:i+1,:), v_f_syn(i,:), t_f(i), t_f(i+1), nro_T, intervals, cr3bp);
end
x_LR_syn(size(x_LR_syn,1) + 1, :) = zeros(1,6);

%Lambert LR (SYN)
 for i=1:np-1
   [Transfer_LR(i).it, it_LR(i), err_LR(i).tr] = Lambert_cycle(x_LR_syn(i:i+1,:), state_time(t_f(i),nro_T), phi0, TOF(i), toll, it_max, cr3bp, options);
 end
toc()

 %% Building Hybrid Solution
 
 output.check_CW = 0;
 output.check_SL = 0;
 output.check_LR = 0;

for i=1:np-1
 if it_CW(i) == it_max
    sprintf('CW Transfer %d failed to converge',i)
     output.check_CW = 1;

     if it_SL(i) == it_max
        sprintf('SL Transfer %d failed to converge',i)
          output.check_SL = 1;
         
           if it_LR(i) == it_max
              sprintf('All methods failed to converge for transfer %d',i)
              output.check_LR = 1;
              return
              %error('First Guesses are not close enough to the Lambert''s solution')
           else
               Lam_syn_hyb(i).transfer = Transfer_LR(i).it(end);
               sprintf('LR Transfer %d converged',i)    
           end
     else
          Lam_syn_hyb(i).transfer = Transfer_SL(i).it(end);
          sprintf('SL Transfer %d converged',i)
    
     end
      
 else
     Lam_syn_hyb(i).transfer = Transfer_CW(i).it(end);
      sprintf('CW Transfer %d converged',i)
 end

end

 
%% Lambert Solution in LVLH Frame 
for i=1:np-1
    Lam_LVLH(i).tv = Lam_syn_hyb(i).transfer.tv;
    Lam_LVLH(i).yv = rsyn2rlvlh_v(Lam_syn_hyb(i).transfer.tv, Lam_syn_hyb(i).transfer.yv(:,1:6),Lam_syn_hyb(i).transfer.yv(:,7:12), cr3bp.mu);
    POSITION_ERROR(i) = norm(Lam_LVLH(i).yv(end,1:3) - r_f(i+1,:));
end

 
%% deltaV calculation

%initial impulse from orbit
dV1 = norm(Lam_syn_hyb(i).transfer.yv(1,4:6) - x0_rel(4:6)');
 
%final brake for each hold point
for i=1:np-1
   dV2(i) =  norm(Lam_syn_hyb(i).transfer.yv(end,4:6));
end
 
%initial impulse from hold point
if np > 2
 for i=2:np-1
    dV3(i-1) =  norm(Lam_syn_hyb(i).transfer.yv(1,4:6));
 end
else
    dV3 = 0;
end
 
dV_tot_adim = dV1 + sum(dV2) + sum(dV3);
dV_tot_dim = dV_tot_adim * cr3bp.L/cr3bp.T*2*pi; %km/s


%% Outputs

%Lambert Hybrid Solution in Synodic
output.Lam_syn = Lam_syn_hyb;

%Lambert Hybrid Solution in LVLH
output.Lambert_LVLH = Lam_LVLH;

%Clohessy-Wiltshire
output.CW.state_syn = x_CW_syn;
output.CW.state_LVLH = x_CW_LVLH;
output.CW.transfer_LVLH = state_CW;
output.CW.Lambert_syn.transfer = Transfer_CW;
output.CW.Lambert_syn.iterations = it_CW;
output.CW.Lambert_syn.error = err_CW;

%Straight Line
output.SL.state_syn = x_SL_syn;
output.SL.state_LVLH = x_SL_LVLH;
output.SL.transfer_LVLH = state_SL;
output.SL.Lambert_syn.transfer = Transfer_SL;
output.SL.Lambert_syn.iterations = it_SL;
output.SL.Lambert_syn.error = err_SL;

%Linearized Relative
output.LR.state_syn = x_LR_syn;
%output.LR.transfer_syn = state_LR;
output.LR.Lambert_syn.transfer = Transfer_LR;
output.LR.Lambert_syn.iterations = it_LR;
output.LR.Lambert_syn.error = err_LR;

%deltaV
output.dV.departure = [dV1 dV3];
output.dV.brake = dV2;
output.dV.total_dim = dV_tot_dim;

 
%% plotting
if plot_yes==1
   plot_1   %linear models
   plot_3   %Lambert
end

 %la domanda da porsi è: utilizzo le eqz adim o dim? Privilegio la rapidità
 %del codice o la precisione numerica? Rispondere dopo. Adesso
 %privilegiare l'omogeneità del codice.
 