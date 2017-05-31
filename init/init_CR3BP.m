function cr3bp = init_CR3BP(m1, m2, params)
% INIT_CR3BP initializes a given CRTBP system.
% 
% CR3BP = INIT_CR3BP(M1, M2, PARAMS) initializes a CRTBP whose primaries
% are defined by their names in string (M1, M2), M1 being the big primary 
% and M2 the small one (e.g. M1 = 'EARTH', M2 = 'MOON').
%
% The fields (all in SI units) computed by this routine include:
%   - CR3BP.mu: the mass ratio [-].
%   - CR3BP.L:  the unit of distance [km].
%   - CR3BP.T:  the unit of time [s].
%
% See code for details.
% 
% BLB 2016

%Primaries init
cr3bp.m1 = init_body(m1);
cr3bp.m2 = init_body(m2);

%Inner parameters
cr3bp.mu =  cr3bp.m2.M/( cr3bp.m1.M + cr3bp.m2.M );         % [-]  µ = m2/(m1 + m2)
cr3bp.L  =  cr3bp.m2.a;                                     % [km] Distance parameter = semi major axis of m2
cr3bp.T  =  cr3bp.m2.T;                                     % [s]  Time parameter = sidereal period of m2
cr3bp.R1 = cr3bp.m1.Req;                                    % [km] m1 radius
cr3bp.R2 = cr3bp.m2.Req;                                    % [km] m2 radius
cr3bp.rh = (cr3bp.mu/3)^1/3;                                % [-]  Hill's radius (adim formula)
cr3bp.libp_precision = params.libp.precision;               % [-]  precision used to compute the libration points Li

%Li initialization
cr3bp.l1 = init_libp(cr3bp, 1, cr3bp.libp_precision);
cr3bp.l2 = init_libp(cr3bp, 2, cr3bp.libp_precision);
cr3bp.l3 = init_libp(cr3bp, 3, cr3bp.libp_precision);
cr3bp.l4 = init_libp(cr3bp, 4, cr3bp.libp_precision);
cr3bp.l5 = init_libp(cr3bp, 5, cr3bp.libp_precision);

%Name
cr3bp.name =  [cr3bp.m1.name, '+', cr3bp.m2.name];

%Distance to manifold
if(strcmp(m1,'EARTH') && strcmp(m2,'MOON'))
    cr3bp.d_man = 50/cr3bp.L;  % 50km for the Earth-Moon system
else
    cr3bp.d_man = 100/cr3bp.L; % 100km ortherwise
end

%Position of the primaries
cr3bp.m1.pos = [-cr3bp.mu 0 0];
cr3bp.m2.pos = [1-cr3bp.mu 0 0];

end

%% Lib points
%--------------------------------------------------------------------------
% Initializes libration point
%--------------------------------------------------------------------------
% @return the structure libp
function libp = init_libp(cr3bp, number, LIBRATION_POINT_PRECISION)
%Number
libp.number = number;

switch number  
    
    case 1
        %gamma_i: distance to the closest primary
        gamma_i = cr3bp.rh - 1/3.0*cr3bp.rh^2- 1/9*cr3bp.rh^3;                    %initial guess
        gamma_i = rtnewt(gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);   %newton-raphson method
        libp.gamma_i = gamma_i;
        libp.c1      = (cr3bp.mu-1.0+libp.gamma_i)/libp.gamma_i;
        
        %Position
        libp.position(1) = 1 - cr3bp.mu - gamma_i;
        libp.position(2) = 0;
        libp.position(3) = 0;
              
    case 2
        %gamma_i: distance to the closest primary
        gamma_i = cr3bp.rh + 1/3.0*cr3bp.rh^2- 1/9*cr3bp.rh^3;
        gamma_i = rtnewt(gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp.gamma_i = gamma_i;
        libp.c1      = (cr3bp.mu-1.0-libp.gamma_i)/libp.gamma_i;
        
        %Position
        libp.position(1) = 1 - cr3bp.mu + gamma_i;
        libp.position(2) = 0;
        libp.position(3) = 0;
                 
    case 3
        %gamma_i: distance to the closest primary
        gamma_i = 7/12.0*cr3bp.mu + 237^2/12^4*cr3bp.mu^3;
        gamma_i = rtnewt(gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp.gamma_i = 1-gamma_i;  %BEWARE  for L3, gamma3 = L3-M1 distance != L3-M2
        libp.c1      = (cr3bp.mu+libp.gamma_i)/libp.gamma_i;
         
        %Position
        libp.position(1) = -1 - cr3bp.mu + gamma_i;
        libp.position(2) = 0;
        libp.position(3) = 0;
                
    case 4
        %gamma_i: distance to the closest primary
        libp.gamma_i = 1;
        
        %Position
        libp.position(1) = -cr3bp.mu + 0.5;
        libp.position(2) = sqrt(3)/2.0;
        libp.position(3) = 0;
               
    case 5
        %gamma_i: distance to the closest primary
        libp.gamma_i = 1;
        
        %Position
        libp.position(1) = -cr3bp.mu + 0.5;
        libp.position(2) = -sqrt(3)/2.0;
        libp.position(3) = 0;        
end

%Energy & Jacobi constant
libp.Ci = jacobi(libp.position, cr3bp.mu);
libp.Ei = -0.5*libp.Ci;
end
 
%-------------------------------------------------------------------------- 
% Using the Newton-Raphson method, find the root of a function known to lie 
% close to y. The root rtnewt will be refined until its accuracy is known 
% within ± precision_gg. polynomialLi is a user-supplied routine that returns both the 
% function value and the first derivative of the function at the point gg_var.
%--------------------------------------------------------------------------
function root = rtnewt(y, precision_gg, mu, number)
gg_var = y;
while(true)
    [f00, df0] = polynomialLi(mu, number, gg_var);
    
    if(abs(f00) < precision_gg);
        break;
    end
    
    gg_var = gg_var - f00/df0;
end
root = gg_var;
end

%--------------------------------------------------------------------------  
% Provides the function value and its first derivative for the 
% newton-raphson method. f corresponds to the equation satisfied by the 
% Li-m2 distance for the L1/L2 cases and by 1-(Li-m1 distance) for the L3 
% case.
%-------------------------------------------------------------------------- 
function [f, df] = polynomialLi(mu, number, y)

switch number
    case 1
        f =  y^5   - (3.0-mu)*y^4 + (3-2*mu)*y^3 - mu*y^2 +  2*mu*y - mu;
        df = 5*y^4 - 4*(3.0-mu)*y^3 + 3*(3-2*mu)*y^2 - 2*mu*y    +  2*mu;
    case 2
        f =  y^5   + (3.0-mu)*y^4 + (3-2*mu)*y^3 - mu*y^2 -  2*mu*y - mu;
        df = 5*y^4 + 4*(3.0-mu)*y^3 + 3*(3-2*mu)*y^2 - 2*mu*y    -  2*mu;
    case 3
        f= y^5 + (7+mu)*y^4 + (19+6*mu)*y^3 -(24+13*mu)*y^2 + (12+14*mu)*y -7*mu;
        df= 5*y^4 + 4*(7+mu)*y^3 + 3*(19+6*mu)*y^2 -2*(24+13*mu)*y + (12+14*mu);
end

end

%% Bodies
%--------------------------------------------------------------------------
% Initializes celestial bodies
%--------------------------------------------------------------------------
function body = init_body(m1)

days = 86400;
body.name = m1;
switch m1
    
    case 'MERCURY'
        
        %Physical parameters
        body.Req = 2439.7;        %[km]
        body.Rm = 2439.7;         %[km]
        body.M = 0.330104e24;     %[kg]
        body.GM = 22032;          %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 57.91e6;         %[kg]
        body.T = 87.9691*days;     %[s]
        
    case 'VENUS'
        
        %Physical parameters
        body.Req = 6051.8;        %[km]
        body.Rm = 6501.8;         %[km]
        body.M = 4.86732e24;      %[kg]
        body.GM = 324858.63;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 108.21e6;        %[km]
        body.T = 224.701*days;     %[s]
        
    case 'EARTH'
        
        %Physical parameters
        body.Req = 6378.14;         %[km]
        body.Rm  = 6371.00;         %[km]
        body.M   = 5.97219e24;      %[kg]
        body.GM  = 398600.440;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 149.60e6;          %[km]
        body.T = 365.25636*days;     %[s]
        
    case 'MOON'
        
        %Physical parameters
        body.Req = 1737.5;                     %[km]
        body.Rm  = 1737.5;                     %[km]
        body.M   = 0.07345814120628661e24;     %[kg] TO BE CONSISTENT WITH HARD CODED VALUE OF mu(Earth-Moon)=0.012150581623434
        body.GM  = 4902.801;                   %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 384400;           %[km]
        body.T = 27.321582*days;    %[s]
        
    case 'MARS'
        
        %Physical parameters
        body.Req = 3396.19;       %[km]
        body.Rm = 3389.50;        %[km]
        body.M = 0.641693e24;     %[kg]
        body.GM = 42828.3;        %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 227.92e6;       %[kg]
        body.T = 686.98*days;     %[s]
        
    case 'JUPITER'
        
        %Physical parameters
        body.Req =  71492;      %[km]
        body.Rm = 69911;        %[km]
        body.M = 1898.13e24;     %[kg]
        body.GM = 126686511;       %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 778.57e6;       %[kg]
        body.T = 4332.589*days;     %[s]
        
    case 'SATURN'
        
        %Physical parameters
        body.Req =  60268;    %[km]
        body.Rm = 58232;       %[km]
        body.M = 568.319e24;     %[kg]
        body.GM = 37931207.8;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 1433.53e6;       %[kg]
        body.T = 10759.22*days;     %[s]
        
    case 'URANUS'
        
        %Physical parameters
        body.Req = 25559;      %[km]
        body.Rm = 25362;        %[km]
        body.M =  86.8103e24;    %[kg]
        body.GM =  5793966;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a =  2872.46e6;      %[kg]
        body.T =  30685.4*days;   %[s]
        
    case 'NEPTUNE'
        
        %Physical parameters
        body.Req = 24764;      %[km]
        body.Rm = 24622;        %[km]
        body.M = 102.410e24;     %[kg]
        body.GM =  6835107;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a =  4495.06e6;      %[kg]
        body.T =  60189*days;    %[s]
        
    case 'PLUTO'
        
        %Physical parameters
        body.Req =  1195;     %[km]
        body.Rm =  1195;       %[km]
        body.M = .01309e24;     %[kg]
        body.GM =  872.4;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a =  5906.38e6;      %[kg]
        body.T =  90465*days;    %[s]
        
    case 'SUN'
        
        %Physical parameters
        body.Req = 696342;                %[km]
        body.Rm =  696342;                %[km]
        body.M  = 1988544e24;             %[kg]
        body.GM = 1.3271244004193938e11;  %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 0;    %[kg]
        body.T = 0;    %[s]
        
    case 'EARTH_AND_MOON'
        %Equivalent mass of the Earth+Moon system based at the center of mass
        %additionnal physical properties are those of the Earth for consistency)
        
        %Physical parameters
        body.Req = 6378.14;        %[km]
        body.Rm = 6371.00;         %[km]
        body.M = 5.97219e24+0.07342e24;       %[kg]
        body.GM = 398600.440+4902.801;      %[km^3.s^-2]
        
        %Orbital parameters
        body.a = 149.60e6;          %[km]
        body.T = 365.25636*days;     %[s]
end
end
