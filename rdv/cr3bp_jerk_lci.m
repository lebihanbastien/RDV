function out = cr3bp_jerk_lci(t,y,mu)
% cr3bp_jerk_lci the jerk of associated to the CR3BP dynamics, in the
% Earth-Moon case. The reference frame for this computation is the LCI, the
% Lunar-Centered Inertial frame.
% 
% BLB 2016

%Output declaration
out = (1:3)';

%--------------------------------------------------------------------------
%We separate y = [r v]
%--------------------------------------------------------------------------
r = y(1:3);
v = y(4:6);

%--------------------------------------------------------------------------
% Position of the primaries in LCI
%--------------------------------------------------------------------------
% Position of the Earth
rlci_1 = -[cos(t) ; sin(t)  ; 0];
% Position of the Moon is [0; 0; 0], so we can get rid of it

%--------------------------------------------------------------------------
% Velocity of the primaries in LCI
%--------------------------------------------------------------------------
% Earth case
vlci_1 = [sin(t) ; -cos(t) ; 0];
% Moon case is still zero!

%--------------------------------------------------------------------------
%Jerk of the moon in barycenter-centered inertial coordinates (BCI)
%--------------------------------------------------------------------------
jbci_2 = +(1-mu)*[sin(t) ; -cos(t)  ; 0];

%--------------------------------------------------------------------------
% Distances between the primaries and the current state y
%--------------------------------------------------------------------------
d1 = norm(rlci_1 - y(1:3));
d2 = norm(y(1:3));

%--------------------------------------------------------------------------
%Phase space derivatives (last three components)
%--------------------------------------------------------------------------
for n = 1:3
    out(n) = + 3*(1-mu)/d1^5 * (r(n) - rlci_1(n))*dot(v - vlci_1, r - rlci_1)...
             -   (1-mu)/d1^3 * (v(n) - vlci_1(n))...
             + 3*(mu)/d2^5   * (r(n))*dot(v, r)...
             -   (mu)/d2^3   * (v(n)) - jbci_2(n);
end
    

end
