function [zrlvlh] = rsyn2rlvlh(t, zrsyn, zrsyn_t, mu)
% Change of coordinates: from RELATIVE Earth-Moon synodic LVLH.
%
% Inputs:
%  - zrsyn a 6*1 vector giving the relative state.
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
% COC: from RSYN to RCLI for the relative state
%--------------------------------------------------------------------------
zrlci = rsyn2rlci(t, zrsyn);

%--------------------------------------------------------------------------
% COC: from RLCI to LVLH for the relative state
%--------------------------------------------------------------------------
zrlvlh = rlci2rlvlh(t, zrlci, zlci_t, mu);

end


