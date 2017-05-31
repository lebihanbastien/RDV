function [zlvlh] = lci2lvlh(t, zlci, zlci_t, mu)
% Change of coordinates: from LCI to LVLH  coordinates.

%--------------------------------------------------------------------------
% Rotation matrix
%--------------------------------------------------------------------------
Cs = lvlh2lciRotMat(zlci_t, t, mu);

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
zlvlh = Cs\(zlci - zlci_t);

end