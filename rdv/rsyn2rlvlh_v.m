function [zlvlh] = rsyn2rlvlh_v(t, zrsyn, zsyn_t, mu)
% Change of coordinates: from RELATIVE Earth-Moon synodic LVLH. Vector
% format.
%
% Inputs:
%  - zrsyn a n*6 matrix giving the relative state.
%  - zrsyn_t a n*6 vector giving the target state. 
%
% Outputs:
%  - zrlvlh a n*6 vector giving the LVLH state.
%
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlvlh = zeros(size(zrsyn));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlvlh(n,:) = rsyn2rlvlh(t(n), zrsyn(n,:)', zsyn_t(n,:)', mu)';
end

end




