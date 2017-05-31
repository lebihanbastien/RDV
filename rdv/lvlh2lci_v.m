function [zlci] = lvlh2lci_v(t, zlvlh, zlci_t, mu)
% Change of coordinates: from LVLH  to LCI  coordinates (vector format).
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlci = zeros(size(zlvlh));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlci(n,:) = lvlh2lci(t(n), zlvlh(n,:)', zlci_t(n,:)', mu)';
end

end




