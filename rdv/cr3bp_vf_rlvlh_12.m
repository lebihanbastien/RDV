function out = cr3bp_vf_rlvlh_12(t,y,mu)
% cr3bp_vf_rlvlh_12 provides the equations of motion for the CR3BP pb seen
% in the LVLH reference frame of a certain target zt. The state y is given
% defined as follows:
% 
%   y(1:6) = the current state z, in the LVLH frame.
%  y(7:12) = the state zt of the target, in the LCI frame, 
%  necessary to define the LVLH frame at each time t.
%
% Contrary to cr3bp_vf_lvlh_12, this routine makes use of the relative
% state in LCI coordinates (hence RLCI coordinates).
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
zlci_t_dot = cr3bp_vf_lci_6(t, zlci_t, mu);
% Store in the output
out(7:12) = zlci_t_dot;

%--------------------------------------------------------------------------
% Computing the necessary matrices for COC, depending on zlci_t.
%--------------------------------------------------------------------------
[Cs, Csdot] = lvlh2lciRotMatdot(zlci_t, t, mu);

%--------------------------------------------------------------------------
% Dynamics ot the current state, in RLCI frame
%--------------------------------------------------------------------------
% Back in RLCI coordinates
zrlci = Cs*y(1:6);

% Dynamics using the LCI vector field
zlcidot = cr3bp_vf_rlci_12(t, [zrlci ; zlci_t], mu);

%--------------------------------------------------------------------------
% Inverse COC to get the vector field in LVLH frame.
%-------------------------------------------------------------------------- 
out(1:6) = Cs\(zlcidot(1:6) - Csdot*y(1:6));

end
