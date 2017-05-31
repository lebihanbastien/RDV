function [Phi_rr, Phi_rv, Phi_vr, Phi_vv] = Phi_RDV(time,n)

Phi_rr = [ 4-3*cos(n*time) 0 0 ; 6*(sin(n*time) - n*time) 1 0 ; 0 0 cos(n*time)];
Phi_rv = [1/n*sin(n*time) 2/n*(1-cos(n*time)) 0 ; 2/n*(cos(n*time)-1) 1/n*(4*sin(n*time)-3*n*time) 0 ; 0 0 1/n*sin(n*time)];
Phi_vr = [3*n*sin(n*time) 0 0 ; 6*n*(cos(n*time)-1) 0 0 ; 0 0 -n*sin(n*time)];
Phi_vv = [cos(n*time) 2*sin(n*time) 0 ; -2*sin(n*time) 4*cos(n*time)-3 0; 0 0 cos(n*time)];

end