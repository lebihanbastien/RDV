function out = cr3bp_vf_rlci_12(t,y,mu)
% cr3bp_vf_rlci_12 provide the equations of motion for the CIRCULAR
% RESTRICTED THREE-BODY PROBLEM (CR3BP) in MATLAB ODE format, using the
% Relative Lunar-Centered Inertial (RLCI) coordinates
%
% OUT = cr3bp_vf_rlci_12(T, Y, MU) computes the first-order
% equations of motion of the CR3BP at time T and state Y, in RLCI
% coordinates.
% The CR3BP mass ratio is MU.
% The state y is given
% defined as follows:
% 
%   y(1:6) = the current state z, in the RLCI frame.
%  y(7:12) = the state zt of the target, in the LCI frame, 
%  necessary to define the LVLH frame at each time t.
%
% Therefore, the absolute state in LCI frame is just y(1:6)+y(7:12)
% 
% BLB 2016

%--------------------------------------------------------------------------
%Output declaration
%--------------------------------------------------------------------------
out = (1:12)';

%--------------------------------------------------------------------------
%Dynamics of the target in the LCI frame (components 7:12)
%--------------------------------------------------------------------------
out(7:12) = cr3bp_vf_lci_6(t, y(7:12), mu);

%--------------------------------------------------------------------------
%Phase space derivatives (components 1:3)
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);

%--------------------------------------------------------------------------
%Phase space derivatives (components 4:6)
%--------------------------------------------------------------------------
% Derivatives in absolute state
temp = cr3bp_vf_lci_6(t, y(1:6)+y(7:12), mu);

% Back to relative state:
out(4) = temp(4) - out(10);
out(5) = temp(5) - out(11);
out(6) = temp(6) - out(12);

end
