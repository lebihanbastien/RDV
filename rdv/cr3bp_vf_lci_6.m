function out = cr3bp_vf_lci_6(t,y,mu)
% cr3bp_vf_lci_6 provide the equations of motion for the CIRCULAR
% RESTRICTED THREE-BODY PROBLEM (CR3BP) in MATLAB ODE format, using the
% Lunar-Centered Inertial (LCI) coordinates
%
% OUT = cr3bp_vf_lci_6(T, Y, MU) computes the first-order
% equations of motion of the CR3BP at time T and state Y, in LCI
% coordinates.
% The CR3BP mass ratio is MU.
% 
% BLB 2016

%Output declaration
out = (1:6)';

%--------------------------------------------------------------------------
%Phase space derivatives (first three components)
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);

%--------------------------------------------------------------------------
%Acceleration of the moon in barycenter-centered inertial coordinates
%(BCI)
%--------------------------------------------------------------------------
abci_2 = -(1-mu)*[cos(t) ; sin(t)  ; 0];

%--------------------------------------------------------------------------
%Position of the primaries in LCI + distance to the current state
%--------------------------------------------------------------------------
% Position of the Earth
rlci_1 = -[cos(t) ; sin(t)  ; 0];
% Position of the Moon is [0; 0; 0], so we can get rid of it
rlci_m = [0 ; 0 ; 0];

%--------------------------------------------------------------------------
% Distances between the primaries and the current state y
%--------------------------------------------------------------------------
d1 = norm(rlci_1 - y(1:3));
d2 = norm(rlci_m - y(1:3));

%--------------------------------------------------------------------------
%Phase space derivatives (last three components)
%--------------------------------------------------------------------------
out(4) = -(1-mu)/d1^3 * (y(1) - rlci_1(1)) - mu/d2^3 * (y(1) - rlci_m(1)) - abci_2(1);
out(5) = -(1-mu)/d1^3 * (y(2) - rlci_1(2)) - mu/d2^3 * (y(2) - rlci_m(2))- abci_2(2);
out(6) = -(1-mu)/d1^3 * (y(3) - rlci_1(3)) - mu/d2^3 * (y(3) - rlci_m(3)) - abci_2(3);

end
