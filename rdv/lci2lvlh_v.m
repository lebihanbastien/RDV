function [zlvlh] = lci2lvlh_v(t, zlci, zlci_t, mu)
% Change of coordinates: from LCI to LVLH  coordinates (vector format).
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlvlh = zeros(size(zlci));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlvlh(n,:) = lci2lvlh(t(n), zlci(n,:)', zlci_t(n,:)', mu)';
end

end




