function [ out ] = linear_cr3bp_derivatives_6(t,y,cr3bp)
% LINEAR_CR3BP_DERIVATIVES_6 provide the equations of motion for the CIRCULAR
% RESTRICTED THREE-BODY PROBLEM (CR3BP) with linear equations in MATLAB ODE format.
%
% OUT = LINEAR_CR3BP_DERIVATIVES_6(T, Y, MU) computes the first-order
% equations of motion of the CR3BPat time T and state Y.
% The CR3BP mass ratio is MU.
%
% The simplified equations (4.1) of motion are available in the article of
% Breakwell & Brown, entitled :
% "The "HALO" family of three-dimensional periodic orbits in the Earth Moon
% resticted 3-body problem" 1979
% 
%
% See also CR3BP_DERIVATIVES_6, CR3BP_DERIVATIVES_42
% 
% AB & LC 2016
% greatly inspired by BLB 2014
mu = cr3bp.mu;
X_lune = cr3bp.m2.pos(1);

%Output declaration
out = (1:6)';

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
% article : other definition of X 
% x_defined = dist(Earth-moon) - X_needed

out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
X = y(1);
Y = y(2);
Z = y(3);
out(4) = 2*y(5) + X - (1-mu)*(X)/((X^2 + Y^2 + Z^2)^(3/2)) - mu*(X-X_lune)/(((X-X_lune)^2 + Y^2 + Z^2)^(3/2));
out(5) = -2 * y(4) + Y - (1-mu)*Y/(((X)^2 + Y^2 + Z^2)^(3/2)) - mu*Y/(((X-X_lune)^2 + Y^2 + Z^2)^(3/2));
out(6) = -(1-mu)*Z/(((X)^2 + Y^2 + Z^2)^(3/2)) - mu*Z/(((X-X_lune)^2 + Y^2 + Z^2)^(3/2));

%out(4) = 2*y(5) + 3* (y(1)-X_lune) + 3/2 * y(3)^2 - mu* (y(1)-X_lune) /(y(3)^3);
%out(5) = -2 * y(4) - mu * y(2)/(y(3)^3);
%out(6) = -y(3) - mu/(y(3)^2);

%out(4) = 3*y(1) - 2*y(5) - 3/2 * y(3)^2 - mu*y(1)/(y(3)^3) + mu/(y(3)^3) - 3;
%out(5) = -mu*y(2)/(y(3)^3) + 2*y(4);
%out(6) = -y(3) - mu/(y(3)^2);

end




