function [zsyn] = lci2syn_v(t, zlci, mu)
% Change of coordinates: from lunar-centered inertial coordinates to
% Earth-Moon synodic (vector format)
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zsyn = zeros(size(zlci));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zsyn(n,:) = lci2syn(t(n), zlci(n,:)', mu)';
end

end




