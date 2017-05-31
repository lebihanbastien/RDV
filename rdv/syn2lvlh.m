function [zlvlh] = syn2lvlh(t, zsyn, zsyn_t, mu)
% Change of coordinates: from absolute Earth-Moon synodic LVLH.
%
% Inputs:
%  - zsyn a 6*1 vector giving the absolute state.
%  - zsyn_t a 6*1 vector giving the target state. 
%
% Outputs:
%  - zlvlh a 6*1 vector giving the LVLH state.
%

%--------------------------------------------------------------------------
% COC: from SYN to LCI for the Target
%--------------------------------------------------------------------------
zlci_t = syn2lci(t, zsyn_t, mu);

%--------------------------------------------------------------------------
% COC: from SYN to CLI for the absolute state
%--------------------------------------------------------------------------
zlci = syn2lci(t, zsyn, mu);

%--------------------------------------------------------------------------
% COC: from LCI to LVLH for the absolute state
%--------------------------------------------------------------------------
zlvlh = lci2lvlh(t, zlci, zlci_t, mu);

end


