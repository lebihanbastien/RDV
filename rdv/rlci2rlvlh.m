function [zlvlh] = rlci2rlvlh(t, zrlci, zlci_t, mu)
% Change of coordinates: from relative LCI to LVLH  coordinates.

%--------------------------------------------------------------------------
% Rotation matrix
%--------------------------------------------------------------------------
Cs = lvlh2lciRotMat(zlci_t, t, mu);

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
zlvlh = Cs\zrlci;
end
