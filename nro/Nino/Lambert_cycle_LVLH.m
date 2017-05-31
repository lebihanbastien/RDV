function [Lam, it, errv] = Lambert_cycle_LVLH(x0_rel_LVLH, x0_T, phi0, TOF, toll, it_max, cr3bp, options)

 err = 1000;
 it = -1;  
 i = 1;
 
 while err > toll && it < it_max
 [Lam(i).tv, Lam(i).yv] = ode113(@(t,y)cr3bp_vf_lvlh_48_alt(t, y, cr3bp.mu),[0 TOF],[x0_rel_LVLH(1,:), x0_T, phi0']', options);

 Phi_Lam = reshape(Lam(i).yv(end,13:48),6,6);
 Phi_rv = Phi_Lam(1:3,4:6); 
 
 r_obt = Lam(i).yv(end,1:3);
 r_desired = x0_rel_LVLH(2,1:3); 
 
 err = norm(r_obt - r_desired);
 errv(i) = err;
 
 dr0 = (r_obt - r_desired)';    
 dV_Lam = inv(Phi_rv) * (dr0);
 x0_rel_LVLH(1,4:6) = x0_rel_LVLH(1,4:6) - dV_Lam';
 
 i=i+1;
 it=it+1;

 end
 