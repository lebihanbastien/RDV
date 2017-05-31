function [zrsyn] = rlvlh2rsyn_v(t, zrlvlh, zsyn_t, mu)
% Change of coordinates: from LVLH to RELATIVE Earth-Moon synodic. Vector
% format.
%
% Inputs:
%  - zrlvlh a n*6 matrix giving the LVLH state.
%  - zrsyn_t a n*6 vector giving the target state, in Synodic EM frame.
%
% Outputs:
%  - zrlvlh a n*6 vector giving the LVLH state.
%
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zrsyn = zeros(size(zrlvlh));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zrsyn(n,:) = rlvlh2rsyn(t(n), zrlvlh(n,:)', zsyn_t(n,:)', mu)';
end

end




