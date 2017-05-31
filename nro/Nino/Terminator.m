function output = Terminator(choice, LVLH, SYN, times_lin, times_cont, settings, nro_T, phi0, n_transf)

r_LVLH = LVLH.r;
v_LVLH = LVLH.v;
r_syn = SYN.r;
v_syn = SYN.v;
r_next = SYN.r_next;
T = times_lin.T;
tt = times_lin.tt;
tt_2 = times_lin.tt_2;
TOF_cont = times_cont.TOF;
int_cont = times_cont.int;

intervals = settings.Luq_int;
toll = settings.toll;
it_max = settings.it_max;
options = settings.options;
cr3bp = settings.cr3bp;   

choice_help = 0;

%% Clohessy Wiltshire (LVLH)
if choice(1)
   [x_CW_LVLH, ~, ~] = Docking(r_LVLH, v_LVLH, tt, tt_2, intervals, nro_T, cr3bp);
    x_CW_syn = rlvlh2rsyn(tt, x_CW_LVLH', state_time(tt,nro_T)', cr3bp.mu)';
    pp = x_CW_LVLH(end,4:6);
      ff = state_time(tt,nro_T)    ;
%   Initial error     
    err_CW = Lambert_cycle_errors(x_CW_syn, r_next, state_time(tt, nro_T), phi0, T, cr3bp, options);
    sol.CW.err1 = err_CW;

    if  ~choice(2)
%       Lambert CW
        [fg, Transfer, it, err] = Lambert_cycle(x_CW_syn, r_next, state_time(tt, nro_T), phi0, T, toll, it_max, cr3bp, options);
        if it < it_max            
           if symmetric_sol(Transfer) 
%             Output generation
              sol.CW.fg = fg;
              sol.CW.Lam = Transfer;
              sol.CW.it = it;
              sol.CW.err = err;
           else
              disp('Error: Symmetric Solution')
              choice_help = 1;
           end           
        else
          disp('MAX it reached')
          choice_help = 1;
        end
    end    
end
         
%% Straight Line   
if choice(2)
   [x_SL_LVLH, ~, ~] = Docking_SL(r_LVLH, v_LVLH, T);
    x_SL_syn         = rlvlh2rsyn(tt, x_SL_LVLH', state_time(tt,nro_T)', cr3bp.mu)';

% Initial error      
    err_SL = Lambert_cycle_errors(x_SL_syn, r_next, state_time(tt,nro_T), phi0, T, cr3bp, options);
    sol.SL.err1 = err_SL;
   
    if ~choice(1)
%       Lambert SL (SYN)
       [fg, Transfer, it, err] = Lambert_cycle(x_SL_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);
        if it < it_max            
           if symmetric_sol(Transfer)  
%             Output generation
              sol.SL.fg = fg;
              sol.SL.Lam = Transfer;
              sol.SL.it = it;
              sol.SL.err = err;
           else
              disp('Error: Symmetric Solution')
              choice_help = 1;
           end           
        else
          disp('MAX it reached')
          choice_help = 1;
        end
    end    
end
       

%% Luquette
if choice(3)
   [x_LR_syn, ~, ~]  = Docking_Luquette(r_syn, v_syn, tt, tt_2, nro_T, intervals, cr3bp);
  
% Initial error
    err_LR = Lambert_cycle_errors(x_LR_syn, r_next, state_time(tt,nro_T), phi0, T, cr3bp, options);
    sol.LR.err1 = err_LR;

    if ~choice(2)      
%      Lambert LR (SYN)
       [fg, Transfer, it, err] = Lambert_cycle(x_LR_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);

        if it < it_max            
           if symmetric_sol(Transfer)  
%             Output generation
              sol.LR.fg = fg;
              sol.LR.Lam = Transfer;
              sol.LR.it = it;
              sol.LR.err = err;
           else
              disp('Error: Symmetric Solution')
              choice_help = 1;
           end           
        else
          disp('MAX it reached')
          choice_help = 1;
        end
    end    
end


%% All three (hp: if the smallest error does not converge, nor the others)

if choice(1) && choice(2) && choice(3)
   err_min = min([err_CW err_SL err_LR]);
   
   if err_min == err_CW
% Lambert CW
      [fg, Transfer, it, err] = Lambert_cycle(x_CW_syn, r_next, state_time(tt, nro_T), phi0, T, toll, it_max, cr3bp, options);
      if it < it_max    
         sprintf('C-W smallest error for transfer %d, converged', n_transf) 
         if symmetric_sol(Transfer)  
%           Output generation
            sol.final.fg = fg;
            sol.final.Lam = Transfer;
            sol.final.it = it;
            sol.final.err = err;
         else
            disp('Error: C-W Symmetric Solution')
            choice_help = 1;
         end
      else
          sprintf('C-W smallest error for transfer %d, failed to converge', n_transf)
          choice_help = 1;
      end
   
      
   elseif err_min == err_SL  
 %        Lambert CW
          [fg, Transfer, it, err] = Lambert_cycle(x_SL_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);
          if it < it_max    
             sprintf('SL smallest error for transfer %d, converged', n_transf) 
             if symmetric_sol(Transfer)  
%               Output generation
                sol.final.fg = fg;
                sol.final.Lam = Transfer;
                sol.final.it = it;
                sol.final.err = err;
             else
                disp('Error: SL Symmetric Solution') 
                choice_help = 1;
             end
         
         else
          sprintf('SL smallest error for transfer %d, failed to converge', n_transf)
          choice_help = 1;
         end
      
   else
 %     Lambert CW
       [fg, Transfer, it, err] = Lambert_cycle(x_LR_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);
        if it < it_max    
           sprintf('LR smallest error for transfer %d, converged', n_transf)
           if symmetric_sol(Transfer)  
%             Output generation
              sol.final.fg = fg;
              sol.final.Lam = Transfer;
              sol.final.it = it;
              sol.final.err = err;
          else
             disp('Error: LR Symmetric Solution')
             choice_help = 1;
           end
             
       else
          sprintf('LR smallest error for transfer %d, failed to converge', n_transf)
          choice_help = 1;
        end
   end
end

%% Continuation

if choice(4) || choice_help 
    if choice_help
        disp('Continuation Algorithm is here to help')
    end
    
   [fg, Transfer, it, err] = Continuation(T, TOF_cont, int_cont, r_syn, v_syn, r_next, tt, nro_T, intervals, phi0, toll, options, it_max, cr3bp); 

%  Output generation 
   sol.cont.fg = fg;
   sol.cont.Lam = Transfer;
   sol.cont.it = it;
   sol.cont.err = err;
end    
%% Final output
sol.help = choice_help;
output = sol;

end

