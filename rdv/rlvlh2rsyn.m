function [zrsyn] = rlvlh2rsyn(t, zlvlh, zrsyn_t, mu)
% Change of coordinates: from LVLH to RELATIVE Earth-Moon synodic.
%
% Inputs:
%  - zlvlh a 6*1 vector giving the LVLH state.
%  - zrsyn_t a 6*1 vector giving the target state. 
%
% Outputs:
%  - zrlvlh a 6*1 vector giving the LVLH state.
%

%--------------------------------------------------------------------------
% COC: from SYN to LCI for the Target
%--------------------------------------------------------------------------
zlci_t = syn2lci(t, zrsyn_t, mu);

%--------------------------------------------------------------------------
% COC: from LVLH to RCLI for the relative state
%--------------------------------------------------------------------------
zrlci = rlvlh2rlci(t, zlvlh, zlci_t, mu);

%--------------------------------------------------------------------------
% COC: from RLCI to RSYN for the relative state
%--------------------------------------------------------------------------
zrsyn = rlci2rsyn(t, zrlci);

end


