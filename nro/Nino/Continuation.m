function [fg, Lam, it, err] = Continuation(TOF, TOF_cont, int_time, r_syn, v_syn, r_next, tt, nro_T, intervals, phi0, toll, options, it_max, cr3bp)

% Initizialing Variables
int_min  = 5;
int_max  = 500;

% times for continuation
TOF_cont     = TOF_cont*2*pi/cr3bp.T;
int_time     = int_time*2*pi/cr3bp.T;
int_time_min = int_min*2*pi/cr3bp.T;
int_time_max = int_max*2*pi/cr3bp.T;
          
% Luquette with low TOF 
[x_LR_syn, ~, ~]  = Docking_Luquette(r_syn, v_syn, tt, tt + TOF_cont, nro_T, intervals, cr3bp);

 
% Lambert with continuation algorithm
TOF_current = TOF_cont;
ind = 1;       
while TOF_current < TOF
      [fg_temp, Lam_temp, it_temp, err_temp] = Lambert_cycle(x_LR_syn, r_next, state_time(tt,nro_T), phi0, TOF_current, toll, it_max, cr3bp, options);
      x_LR_syn = Lam_temp(1,1:6);
      int_time = min(int_time_max, max(int_time*3/it_temp,int_time_min));
      TOF_current = min(TOF_current + int_time, TOF);
end

fg = fg_temp;    
Lam = Lam_temp;
it = it_temp;
err = err_temp;
    
    
    
    
    
    
    
    
    
%       if  ind >= 3
%           dV_cc = V_cont(ind-1,:) + (TT_dV(ind) - TT_dV(ind-1))/(TT_dV(ind-1) - TT_dV(ind-2))*(V_cont(ind-1,:) - V_cont(ind-2,:));
%           x_LR_syn_cont(i,1:6) = dV_cc;
%       end

%       
%       V_cont(ind,:) = x_LR_syn_cont(i,4:6);
%       TT_dV(ind) = TOF_current;
%       TT_dV(ind+1) = TOF_current;
%       ind = ind+1;
% end
% x_LR_syn_cont(i+1, :) = [r_f_syn(i+1,:), v_f_syn(i+1,:)];
%     
