function [zlci] = rsyn2rlci_v(t, zsyn)
% Change of coordinates: from relative Earth-Moon synodic to relative
% lunar-centered inertial coordinates (vector format)
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlci = zeros(size(zsyn));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlci(n,:) = rsyn2rlci(t(n), zsyn(n,:)')';
end

end




