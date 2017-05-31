function [x0_corr, deltaV, cc] = Docking(r_hold, v_hold, t_i, t_f, intervals, nro_T, cr3bp)

% divide TOF in sub-intervals
times = (t_f - t_i)/intervals;
t_int = t_i:times:t_f;

% target state for distances
x0_T = reshape(state_time(t_int(1:end-1),nro_T),[],6);

% variables
Phi_temp = zeros(6,6,length(t_int)-1);

for i = 1:length(t_int)-1
    %x0_TT(i,:) = state_time(t_int(i),nro_T)

    r2_T = [x0_T(i,1) - 1 + cr3bp.mu, x0_T(i,2), x0_T(i,3)]; % distance from Moon
    a_T = norm(r2_T);
    n = sqrt(cr3bp.mu/(a_T^3));   % "Mean" orbital velocity
    A_CW = [0     0     0     1    0    0
            0     0     0     0    1    0
            0     0     0     0    0    1
            0     0     0     0    0   2*n
            0   -n^2    0     0    0    0
            0     0   3*n^2  -2*n  0    0];  
% STM for each sub-interval
    Phi_temp(:,:,i) = expm(A_CW*times);
end

% build the STM for the whole transfer
Phi = Phi_temp(:,:,end);

for i=1:length(t_int)-2
    Phi = Phi*Phi_temp(:,:,end-i);
end
 
%verify STM is well-conditioned
cc = rcond(Phi);

% transfer 
%[dv0i_p, deltaV] = transfer_rel (A_CW, r_hold(2,:)', TOF, [r_hold(1,:) v_hold]');
[dv0i_p, deltaV] = transfer_rel_Luquette (Phi, r_hold(2,:)', [r_hold(1,:) v_hold]');
x0_corr = [r_hold(1,:) dv0i_p'];

% options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% [~,state] = ode45( @Clohessy_Wiltshire, [0, TOF] ,x0_corr',options,n);

end