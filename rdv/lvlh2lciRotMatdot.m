function [Cs, Csdot, e1, e2, e3, e1d, e2d, e3d, e1dd, e2dd, e3dd] = lvlh2lciRotMatdot(zlci_t, t, mu)
% Rotation matrix for the change of coordinates: from LVLH to LCI coordinates.
% The derivative of Cs is also computed

%--------------------------------------------------------------------------
% Decomposing position/velocity
%--------------------------------------------------------------------------
r = zlci_t(1:3);
v = zlci_t(4:6);

%--------------------------------------------------------------------------
% a = dv/dt in LCI is obtained by applying cr3bp_vf_lci_6
%--------------------------------------------------------------------------
zlci_t_dot = cr3bp_vf_lci_6(t, zlci_t, mu);
a = zlci_t_dot(4:6);

%--------------------------------------------------------------------------
% j = da/dt in LCI is obtained by applying cr3bp_jerk_lci
%--------------------------------------------------------------------------
j = cr3bp_jerk_lci(t,zlci_t,mu);

%--------------------------------------------------------------------------
% Preliminaries: 
% 1. We compute crv = cross(r,v) and cra = cross(r, a), etc
% 2. we compute k = ||r|| and h = ||r x v|| and their derivatives
%--------------------------------------------------------------------------
crv = cross(r, v);
cra = cross(r, a);
cva = cross(v, a);
crj = cross(r, j);

k   = norm(r);
kd  = r'*v/k;
kdd = (v'*v + r'*a - kd^2)/k;


h   = norm(crv);
hd  = crv'*cra/h;
hdd = ( cra'*cra + crv'*(cva + crj) - hd^2)/h;

%--------------------------------------------------------------------------
% Building C = [e1 e2 e3]
%--------------------------------------------------------------------------
% z-vector e3 is -r/||r||
e3 = - r/k;

% y-vector e2 is - r x v/||r x v||
e2 = - crv/h;

% x-vector e1 is e2 x e3
e1 = cross(e2, e3);

% C = [e1 e2 e3]
C = [e1 e2 e3];

%--------------------------------------------------------------------------
% Building Cdot = [de1/dt de2/dt de3/dt]
%--------------------------------------------------------------------------
e3d = -(k*v - kd*r)/k^2;

e2d = -(h*cra - hd*crv)/h^2;

e1d = cross(e2d, e3) + cross(e2, e3d);

% Compute: Cdot = [e1d e2d e3d]
Cdot = [e1d e2d e3d];

%--------------------------------------------------------------------------
% Building Cddot = [d^2 e1/dt^2 d^2 e2/dt^2  d^2 e3/dt^2]
%--------------------------------------------------------------------------
e3dd = -( (2*kd^2 - k*kdd)*r - 2*k*kd*v + k^2*a)/k^3;

e2dd = -( h^2*(cva + crj) - 2*h*hd*cra + (2*hd^2 - h*hdd)*crv)/h^3;

e1dd = cross(e2dd, e3) + 2*cross(e2d, e3d) + cross(e2, e3dd);

% Compute: Cddot = [e1dd e2dd e3dd]
Cddot = [e1dd e2dd e3dd];

%--------------------------------------------------------------------------
% Building Cs = [ C    0 ;
%                  Cdot C];
%--------------------------------------------------------------------------
Cs = [ C    zeros(3)  ;
       Cdot       C  ];
   
%--------------------------------------------------------------------------
% Building Csdot = [ Cdot    0 ;
%                    Cddot Cdot];
%--------------------------------------------------------------------------
Csdot = [ Cdot    zeros(3) ;
          Cddot     Cdot  ];
end



