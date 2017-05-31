function out = cr3bp_derivatives_6(t,y,mu)
% CR3BP_DERIVATIVES_6 provide the equations of motion for the CIRCULAR
% RESTRICTED THREE-BODY PROBLEM (CR3BP) in MATLAB ODE format.
%
% OUT = CR3BP_DERIVATIVES_6(T, Y, MU) computes the first-order
% equations of motion of the CR3BPat time T and state Y.
% The CR3BP mass ratio is MU.
%
% The equations of motion are available in chapter 2 of Koon et al.
% "Dynamical Systems, the Three-Body Problem and Space Mission Design" 2006
% <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% See also CR3BP_DERIVATIVES_42
% 
% BLB 2014

%Output declaration
out = (1:6)';

%--------------------------------------------------------------------------
% Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
%--------------------------------------------------------------------------
d1_ub = d1_u_barre(mu,y(1),y(2),y(3));

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = -d1_ub(1) + 2*y(5);
out(5) = -d1_ub(2) - 2*y(4);
out(6) = -d1_ub(3);

end

%-------------------------------------------------------------------------%
% First derivative of the energy potential
%-------------------------------------------------------------------------%
function xout = d1_u_barre(mu,x,y,z)
xout(1) = (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2)) - ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) - x;
xout(2) = (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - y - (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
xout(3) = (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
end

