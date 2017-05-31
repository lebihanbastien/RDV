function [dv_i_p, deltaV_tot] = transfer_rel_Luquette (Phi, r_hold, initial_state)

Phi_rr = Phi(1:3,1:3);
Phi_rv = Phi(1:3,4:6);    
Phi_vr = Phi(4:6,1:3);
Phi_vv = Phi(4:6,4:6);

dr0 = initial_state(1:3);
dv_i_m = initial_state(4:6);

dv_i_p = Phi_rv \ (r_hold - Phi_rr* dr0);
dv_f_m = (Phi_vr - Phi_vv * (Phi_rv \ Phi_rr)) * dr0 + Phi_vv * (Phi_rv \ r_hold);
dv_f_p = 0;


deltaV_i = norm(dv_i_p - dv_i_m);
deltaV_f = norm(dv_f_p - dv_f_m);

deltaV_tot = abs(deltaV_i) + abs(deltaV_f);