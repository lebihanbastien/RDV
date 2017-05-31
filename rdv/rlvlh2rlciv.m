function [zrlci] = rlvlh2rlciv(t, zlvlh, zlci_t, mu)
% Change of coordinates: from LVLH  to RLCI  coordinates (vector format).
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zrlci = zeros(size(zlvlh));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zrlci(n,:) = rlvlh2rlci(t(n), zlvlh(n,:)', zlci_t(n,:)', mu)';
end

end




