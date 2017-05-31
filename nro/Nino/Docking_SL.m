function [x0_corr, deltaV, state] = Docking_SL(r_hold, v_hold, TOF)

% transfer
A_SL = [0     0     0     1    0    0
        0     0     0     0    1    0
        0     0     0     0    0    1
        0     0     0     0    0    0
        0     0     0     0    0    0
        0     0     0     0    0    0];
    
[dv0i_p, deltaV] = transfer_rel (A_SL, r_hold(2,:)', TOF, [r_hold(1,:) v_hold]');
x0_corr = [r_hold(1,:) dv0i_p'];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~,state] = ode45( @Straight_line, [0, TOF] , x0_corr',options);

end