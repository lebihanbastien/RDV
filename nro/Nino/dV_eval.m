function output = dV_eval(Lam, x0_rel, np, cr3bp)

% Initialization
dV2 = zeros(np-1,1);
dV3 = zeros(np-2,1);  %np-2 because before the first impulse from the orbit, the chaser has a non-zero velocity

% initial impulse from orbit
dV1 = norm(Lam{1}(1,4:6) - x0_rel(4:6)');
 
% final brake for each hold point
for i=1:np-1
    dV2(i) =  norm(Lam{i}(end,4:6));
end
 
% initial impulse from hold point
if np > 2
   for i=2:np-1
       dV3(i-1) =  norm(Lam{i}(1,4:6));
   end
else
    dV3 = 0;
end
 
dV_tot_adim = dV1 + sum(dV2) + sum(dV3);
dV_tot_dim = dV_tot_adim * cr3bp.L/cr3bp.T*2*pi; %km/s

% deltaV
output.dV.departure = [dV1 dV3'];
output.dV.brake = dV2;
output.dV.total_dim = dV_tot_dim;