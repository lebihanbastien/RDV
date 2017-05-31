function [zlci] = syn2lci_v(t, zsyn, mu)
% Change of coordinates: from Earth-Moon synodic to lunar-centered inertial
% coordinates (vector format)
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlci = zeros(size(zsyn));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlci(n,:) = syn2lci(t(n), zsyn(n,:)', mu)';
end

end




