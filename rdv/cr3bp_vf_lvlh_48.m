function out = cr3bp_vf_lvlh_48(t, y, mu)
% cr3bp_vf_lvlh_48 provides the equations of motion for the CR3BP pb seen
% in the LVLH reference frame of a certain target zt. The variational
% equations are included. The state y is defined as follows:
% 
%   y(1:6)   = the current state z, in the LVLH frame.
%  y(7:12)   = the state zt of the target, in the LCI frame, 
%  necessary to define the LVLH frame at each time t.
%   y(13:48) =  the STM in LVLH frame.
% 
% BLB 2016

%--------------------------------------------------------------------------
% Output declaration
%--------------------------------------------------------------------------
out = (1:12)';

%--------------------------------------------------------------------------
% Dynamics of the target in the LCI frame (last 6 components)
%--------------------------------------------------------------------------
% State of the target, in LCI coordinates.
zlci_t = y(7:12);
% Derivatives of dzt/dt, in LCI coordinates.
out(7:12) = cr3bp_vf_lci_6(t, zlci_t, mu);

%--------------------------------------------------------------------------
%Acceleration of the moon in barycenter-centered inertial coordinates (BCI)
%--------------------------------------------------------------------------
ct = cos(t);
st = sin(t);
abci_2 = -(1-mu)*[ct ; st  ; 0];

%--------------------------------------------------------------------------
% First, for the relative state, we can actually force the first 
% three components of out
%-------------------------------------------------------------------------- 
out(1:3) = y(4:6);

%--------------------------------------------------------------------------
% Computing the necessary vectors
%--------------------------------------------------------------------------
[~, ~, e1, e2, e3, e1d, e2d, e3d, e1dd] = lvlh2lciRotMatdot(zlci_t, t, mu);

%--------------------------------------------------------------------------
% Positions of the primary in Rlvlh coordinates
%--------------------------------------------------------------------------
% Positions in LCI coordinates: 
rlci_1 = -[ct ; st ; 0];
rlci_2 =  [0  ; 0  ; 0];

% Back in LVLH
rlvlh_1 = [e1 e2 e3]\(rlci_1 - zlci_t(1:3));
rlvlh_2 = [e1 e2 e3]\(rlci_2 - zlci_t(1:3));

%--------------------------------------------------------------------------
% Distances between the primaries and the current state y
%--------------------------------------------------------------------------
d1 = norm(rlvlh_1 - y(1:3));
d2 = norm(rlvlh_2 - y(1:3));

%--------------------------------------------------------------------------
% Computing the coefficicents of the vector field (x'*y = dot(x,y))
%
% Note that out(10:12)+abci_2 is the acceleration of the Target in Inertial
% coordinates centered at the barycenter of the Earth-Moon system.
%--------------------------------------------------------------------------
b1  = - e1'  * (out(10:12)+abci_2);
b2  = - e2'  * (out(10:12)+abci_2);
b3  = - e3'  * (out(10:12)+abci_2);
b4  = 2*e1d' * e2;
b5  = 2*e1d' * e3;
b6  =   e1d' * e1d;
b7  =   e1dd'* e2;
b8  =   e1dd'* e3;
b9  =   e2d' * e2d;
b10 =   e2d' * e3d;
b11 =   e3d' * e3d;

%--------------------------------------------------------------------------
% Vector field out(4:)
%--------------------------------------------------------------------------
out(4) = b1           + b4*y(5) + b5*y(6) + b6*y(1) + b7*y(2) + b8*y(3)...
         - (1-mu)/d1^3*(y(1) - rlvlh_1(1)) - mu/d2^3*(y(1) - rlvlh_2(1));
out(5) = b2 - b4*y(4)                     - b7*y(1) + b9*y(2) + b10*y(3)...
         - (1-mu)/d1^3*(y(2) - rlvlh_1(2)) - mu/d2^3*(y(2) - rlvlh_2(2));
     
out(6) = b3 - b5*y(4)                     - b8*y(1) + b10*y(2) + b11*y(3)...
         - (1-mu)/d1^3*(y(3) - rlvlh_1(3)) - mu/d2^3*(y(3) - rlvlh_2(3));
     
%--------------------------------------------------------------------------
% STM derivatives: out(13:48)
%--------------------------------------------------------------------------
STM = reshape(y(13:48), [6, 6]);

% Variational equations matrix
A = [ 0    0    0     1    0    0  ;
      0    0    0     0    1    0  ;
      0    0    0     0    0    1  ;
      b6   b7   b8    0    b4   b5 ;
     -b7   b9   b10  -b4   0    0  ;
     -b8   b10  b11  -b5   0    0  ];

%Derivative
dSTM = A * STM;

%Storage in output
out(13:48) = reshape(dSTM, [36,1]);
end
