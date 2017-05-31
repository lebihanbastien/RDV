function [yv_li, T_li] =  third_order_orbit(orbit, t, cst)
% THIRD_ORDER_ORBIT computes a third-order approximation of initial
% conditions position+velocity, for three types of orbits: halo, planar
% lyapunov and vertical lyapunov.
%
% [YV_LI, T_LI] =  THIRD_ORDER_ORBIT(ORBIT, T, CST) computes a third-order
% approximation of initial conditions position+velocity from information
% contained in the structure ORBIT. The initial conditions are given at
% time T, and are stored in YV_LI. Moreover, an estimate of the period of
% the orbit is given in T_LI.
%
% see: "Analytic construction of periodic orbits about the collinear points"
% <a href="matlab: 
% web('http://adsabs.harvard.edu/full/1980CeMec..22..241R','-browser')">(link)</a>
% Richardson 1980
%
% The orbit must have been properly initialized via the routine
% INIT_ORBIT.
%
% See also INIT_ORBIT, RICHARDSON_COEFFICIENTS
%
% BLB 2015

%--------------------------------------------------------------------------
% Init
%--------------------------------------------------------------------------
%Orbit parameters
gamma_i = orbit.li.gamma_i;
mu      = orbit.cr3bp.mu;
dm      = orbit.dm;

%Richardson coefficients
RC = richardson_coefficients(orbit.cr3bp, orbit.li.number);

%Phases
Phi = 0.0; %arbitrary

%X and Z-amplitude in Li-frame
switch(orbit.type)
    
    case cst.orbit.type.HALO
        Az_li = orbit.Az_estimate/gamma_i;
        Ax_li = sqrt(-(RC.l2*Az_li^2 + RC.Delta)/RC.l1);
        
    case cst.orbit.type.PLYAP
        Az_li = 0.0;
        Ax_li = orbit.Ax_estimate/gamma_i;
        
    case cst.orbit.type.VLYAP
        Az_li = orbit.Az_estimate/gamma_i;
        Ax_li = 0.0;
end


%nu parameter
nu = 1+RC.s1*Ax_li^2+RC.s2*Az_li^2;

%--------------------------------------------------------------------------
% Period estimate
%--------------------------------------------------------------------------
T_li = 2*pi/(RC.omega_p*nu);

%--------------------------------------------------------------------------
% State estimate @ t  in the Li-frame
%--------------------------------------------------------------------------
tau1 = RC.omega_p*nu*t + Phi;
dtau1dt = RC.omega_p*nu;

yv_li = richardson_init_cond(Ax_li, Az_li, tau1, dtau1dt, dm, RC);

%--------------------------------------------------------------------------
% Back into CR3BP frame
%--------------------------------------------------------------------------
% Special case of the x dimension
switch(orbit.li.number)
    case 1
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 - gamma_i;
        
    case 2
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 + gamma_i;
        
    case 3
        yv_li(1) = yv_li(1)*gamma_i - mu - gamma_i;
        
    otherwise %L2
        yv_li(1) = yv_li(1)*gamma_i - mu + 1 + gamma_i;
end

% Other dimensions
for i = 2:6
    yv_li(i) = gamma_i*yv_li(i);
end


end


function yv_li = richardson_init_cond(Ax_li, Az_li, tau1, dtau1dt, dm, RC)
% Initial conditions at time TAU1

yv_li(1) = RC.a21*Ax_li^2 + RC.a22*Az_li^2 - Ax_li*cos(tau1) + (RC.a23*Ax_li^2 - RC.a24*Az_li^2)*cos(2*tau1) + (RC.a31*Ax_li^3 - RC.a32*Ax_li*Az_li^2)*cos(3*tau1);
yv_li(2) = RC.kappa*Ax_li*sin(tau1) + (RC.b21*Ax_li^2 - RC.b22*Az_li^2)*sin(2*tau1) + (RC.b31*Ax_li^3 - RC.b32*Ax_li*Az_li^2)*sin(3*tau1);
yv_li(3) = dm*Az_li*cos(tau1) + dm*RC.d21*Ax_li*Az_li*(cos(2*tau1)-3) + dm*(RC.d32*Az_li*Ax_li^2 - RC.d31*Az_li^3)*cos(3*tau1);

yv_li(4) = dtau1dt*( Ax_li*sin(tau1) - 2*(RC.a23*Ax_li^2 - RC.a24*Az_li^2)*sin(2*tau1) -3*(RC.a31*Ax_li^3 - RC.a32*Ax_li*Az_li^2)*sin(3*tau1)  );
yv_li(5) = dtau1dt*( RC.kappa*Ax_li*cos(tau1) + 2*(RC.b21*Ax_li^2 - RC.b22*Az_li^2)*cos(2*tau1) + 3*(RC.b31*Ax_li^3 - RC.b32*Ax_li*Az_li^2)*cos(3*tau1)   );
yv_li(6) = dtau1dt*( -dm*Az_li*sin(tau1) -2*dm*RC.d21*Ax_li*Az_li*sin(2*tau1) - 3*dm*(RC.d32*Az_li*Ax_li^2 - RC.d31*Az_li^3)*sin(3*tau1)  );

end