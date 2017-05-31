function output = Terminator(choice, LVLH, SYN, times_lin, times_cont, settings, nro_T, phi0)

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

%% Clohessy Wiltshire (LVLH)
if choice(1)
   [x_CW_LVLH, ~, cc_CW] = Docking(r_LVLH, v_LVLH, tt, tt_2, intervals, nro_T, cr3bp);
    x_CW_syn = rlvlh2rsyn(tt, x_CW_LVLH', state_time(tt,nro_T)', cr3bp.mu)';
        
% Lambert CW
   [fg, Transfer, it, err] = Lambert_cycle(x_CW_syn, r_next, state_time(tt, nro_T), phi0, T, toll, it_max, cr3bp, options);
   
% for error     
%    err_CW = Lambert_cycle_errors(x_CW_syn, r_next, state_time(tt, nro_T), phi0, T, cr3bp, options);

%    sol.CW.err = err_CW;
%    sol.CW.cc = cc_CW;
% Output generation
    sol.CW.fg = fg;
    sol.CW.Lam = Transfer;
    sol.CW.it = it;
    sol.CW.err = err;
end
         
         
%% Straight Line   
if choice(2)
  [x_SL_LVLH, ~, ~] = Docking_SL(r_LVLH, v_LVLH, T);
   x_SL_syn         = rlvlh2rsyn(tt, x_SL_LVLH', state_time(tt,nro_T)', cr3bp.mu)';
        
% Lambert SL (SYN)
  [fg, Transfer, it, err] = Lambert_cycle(x_SL_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);
%    sol.SL.err = Lambert_cycle_errors(x_SL_syn, r_next, state_time(tt,nro_T), phi0, T, cr3bp, options);


% Output generation
    sol.SL.fg = fg;
    sol.SL.Lam = Transfer;
    sol.SL.it = it;
    sol.SL.err = err;
end
       

%% Luquette
if choice(3)      
  [x_LR_syn, ~, cc_LR]  = Docking_Luquette(r_syn, v_syn, tt, tt_2, nro_T, intervals, cr3bp);
  [fg, Transfer, it, err] = Lambert_cycle(x_LR_syn, r_next, state_time(tt,nro_T), phi0, T, toll, it_max, cr3bp, options);
  
% for error  
%    err_LR = Lambert_cycle_errors(x_LR_syn, r_next, state_time(tt,nro_T), phi0, T, cr3bp, options);

%     sol.LR.err = err_LR;
%     sol.LR.cc = cc_LR;


% Output generation 
    sol.LR.fg = fg;
    sol.LR.Lam = Transfer;
    sol.LR.it = it;
    sol.LR.err = err;
end


%% Continuation

if choice(4)
   [fg, Transfer, it, err] = Continuation(T, TOF_cont, int_cont, r_syn, v_syn, r_next, tt, nro_T, intervals, phi0, toll, options, it_max, cr3bp); 

% Output generation 
   sol.cont.fg = fg;
   sol.cont.Lam = Transfer;
   sol.cont.it = it;
   sol.cont.err = err;
end    
%% Final output
output = sol;

end

