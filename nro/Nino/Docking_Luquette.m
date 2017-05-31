function [x0_corr, deltaV, cc] = Docking_Luquette(r_hold, v_hold, t_i, t_f, nro_T, intervals, cr3bp)

% divide TOF in sub-intervals
times = (t_f - t_i)/intervals;
t_int = t_i:times:t_f;

% target state for distances
x0_T = reshape(state_time(t_int(1:end-1),nro_T),[],6);

% costants
mu_1 = (1 - cr3bp.mu);
mu_2 = cr3bp.mu;
n1 = [0 -1 0; 1 0 0; 0 0 0 ];
nn = [-1 0 0; 0 -1 0; 0 0 0];

% variables
Phi_temp = zeros(6,6,length(t_int)-1);

for i=1:length(t_int)-1
% distances Target from Moon and Earth    
    r1 = [x0_T(i,1) + cr3bp.mu, x0_T(i,2), x0_T(i,3) ];  %Earth
    r2 = [x0_T(i,1) - 1 + cr3bp.mu, x0_T(i,2), x0_T(i,3)]; %Moon
    r1_T = norm(r1);
    r2_T = norm(r2);
% Luquette
    csi = - (mu_1/r1_T^3 + mu_2/r2_T^3)*eye(3) + 3 * mu_1/r1_T^5 * (r1*r1') + 3*mu_2/r2_T^5 * (r2*r2');
    A_LRT = [zeros(3,3),  eye(3); csi - nn, -2*n1];
% STM for each sub-interval
    Phi_temp(:,:,i) = expm(A_LRT*times);

end

% build the STM for the whole transfer
Phi = Phi_temp(:,:,end);

for i=1:length(t_int)-2
    Phi = Phi*Phi_temp(:,:,end-i);
end

%verify STM is well-conditioned
cc = rcond(Phi);

%transfer
[dv0i_p, deltaV] = transfer_rel_Luquette (Phi, r_hold(2,:)', [r_hold(1,:) v_hold]');
x0_corr = [r_hold(1,:) dv0i_p'];

% A_LRT_vec = reshape(A_LRT,[],1);
% options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% [t_int,state] = ode45( @LRT, [0, t_f-t_i], x0_corr', options, A_LRT_vec);