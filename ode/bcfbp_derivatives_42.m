function out = bcfbp_derivatives_42(t,y,mu, omega0)
% BCFBP_DERIVATIVES_42 provide the equations of motion for the EARTH-MOON 
% BICIRCULAR FOUR-BODY PROBLEM (BCP) in MATLAB ODE format.
%
% OUT = BCFBP_DERIVATIVES_42(T, Y, MU, OMEGA0) computes the first-order
% equations of motion of the EARTH-MOON BCP at time T and state Y.
% The Earth-Moon mass ratio is MU. The initial position of the fourth body 
% (e.g. the Sun) is given by the scalar OMEGA0. On top of the
% 6-dimensionnal state, the 36 equations that govern the evoluation of the
% State Transition Matrix (STM) are also provided, for a total of 42
% first-order equations of motion.
%
% The equations of motion are available in equation (5.3.1) of Koon et al.
% "Dynamical Systems, the Three-Body Problem and Space Mission Design" 2006
% <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% Note: the Sun parameters (mass, mean motion, semi-major axis) are hard
% coded in the routine.
%
% TODO: Outsource the hard coded values of ms, as, omegaS
%
% See also BCFBP_DERIVATIVES_6 
%
% BLB 2015

ms = 328900.54; 
as = 388.81114;
omegaS = -0.925195985520347;

%Output declaration
out = (1:42)';

%--------------------------------------------------------------------------
% Update first derivatives of the potential \bar{U} (cf Koon et al. 2006)
%--------------------------------------------------------------------------
d1_ub = d1_u_b4(mu, ms, as, y(1),y(2),y(3), omega0+omegaS*t);
d2_ub = d2_u_b4(mu, ms, as, y(1),y(2),y(3), omega0+omegaS*t);

%--------------------------------------------------------------------------
%Phase space derivatives
%--------------------------------------------------------------------------
out(1) = y(4);
out(2) = y(5);
out(3) = y(6);
out(4) = d1_ub(1) + 2*y(5);
out(5) = d1_ub(2) - 2*y(4);
out(6) = d1_ub(3);

%--------------------------------------------------------------------------
% STM derivatives
%--------------------------------------------------------------------------
%STM is updated
Phi = eye(6);
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        Phi(i,j) = y(6+m);
    end
end

% Omega 3*3
Omega = zeros(3);
Omega(1,2) = 1;
Omega(2,1) = -1;

%Matrices triviales
I3 = eye(3);
Z3 = zeros(3);

%Super matrix Df(x) (6*6)
Df = [Z3 I3 ; d2_ub 2*Omega];

%Derivative
d_Phi = Df * Phi;

%Storage in output
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        out(6+m) = d_Phi(i,j);
    end
end



end

%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function xout = d1_u_b4(mu, ms, as, x,y,z, theta)
% First derivatives of the energy potential
%
% BLB 2015

r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
rs = sqrt((x-as*cos(theta))^2 + (y-as*sin(theta))^2 + z^2);
xout(1) = x - (1-mu)/r1^3*(x+mu) - mu/r2^3*(x-1+mu) - ms/rs^3*(x-as*cos(theta)) - ms/as^2*cos(theta);
xout(2) = y - (1-mu)/r1^3*y      - mu/r2^3*y        - ms/rs^3*(y-as*sin(theta)) - ms/as^2*sin(theta);
xout(3) =   - (1-mu)/r1^3*z      - mu/r2^3*z        - ms/rs^3*z;
end


function xout = d2_u_b4(mu, ms, as, x,y,z, theta)
% Seconde derivatives of the energy potential
%
% BLB 2015

xout(1,1) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - ms/((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*ms*(2*x - 2*as*cos(theta))^2)/(4*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
xout(1,2) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + (3*ms*(2*x - 2*as*cos(theta))*(2*y - 2*as*sin(theta)))/(4*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2));
xout(1,3) = (3*ms*z*(2*x - 2*as*cos(theta)))/(2*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) + (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));

xout(2,1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + (3*ms*(2*x - 2*as*cos(theta))*(2*y - 2*as*sin(theta)))/(4*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2));
xout(2,2) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - ms/((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*ms*(2*y - 2*as*sin(theta))^2)/(4*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
xout(2,3) = (3*ms*z*(2*y - 2*as*sin(theta)))/(2*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);

xout(3,1) = (3*ms*z*(2*x - 2*as*cos(theta)))/(2*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) + (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
xout(3,2) = (3*ms*z*(2*y - 2*as*sin(theta)))/(2*((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2)) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);
xout(3,3) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - ms/((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + (3*ms*z^2)/((x - as*cos(theta))^2 + z^2 + (y - as*sin(theta))^2)^(5/2);

end


