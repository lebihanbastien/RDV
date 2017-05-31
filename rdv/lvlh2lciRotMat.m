function [Cs] = lvlh2lciRotMat(zlci_t, t, mu)
% Rotation matrix for the change of coordinates: from LVLH to LCI coordinates.

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
% Preliminaries: 
% 1. We compute crv = cross(r,v) and cra = cross(r, a)
% 2. we compute k = ||r|| and h = ||r x v|| and their derivatives
%--------------------------------------------------------------------------
crv = cross(r, v);
cra = cross(r, a);

k  = norm(r);
kd = dot(r, v)/k;

h  = norm(crv);
hd = dot(crv, cra)/h;


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

% Cdot = [e1d e2d e3d]
Cdot = [e1d e2d e3d];

%--------------------------------------------------------------------------
% Building Cs = [ C    0 ;
%                  Cdot C];
%--------------------------------------------------------------------------
Cs = [ C    zeros(3) ;
       Cdot       C  ];
end



