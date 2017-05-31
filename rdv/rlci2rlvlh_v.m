function [zlvlh] = rlci2rlvlh_v(t, zrlci, zlci_t, mu)
% Change of coordinates: from Relative LCI to LVLH  coordinates 
% (vector format).
%--------------------------------------------------------------------------
% Prealloc
%--------------------------------------------------------------------------
zlvlh = zeros(size(zrlci));

%--------------------------------------------------------------------------
% COC
%--------------------------------------------------------------------------
for n = 1:size(t,1)   
    zlvlh(n,:) = rlci2rlvlh(t(n), zrlci(n,:)', zlci_t(n,:)', mu)';
end

end




