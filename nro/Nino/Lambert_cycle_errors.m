function err = Lambert_cycle_errors(x0_rel, r_next, x0_T, phi0, TOF, cr3bp, options)

[~, Lam] = ode113(@(t,y)crtbp(t,y,cr3bp.mu),[0 TOF],[x0_rel, x0_T, phi0']', options); 
 
r_obt = Lam(end,1:3);
r_desired = r_next; 
err = norm(r_obt - r_desired);
 





