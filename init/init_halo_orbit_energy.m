function orbit = init_halo_orbit_energy(cr3bp, li, family, Cjac, cst)
% INIT_HALO_ORBIT_ENERGY(CR3BP, LI, FAMILY, CJAC, CST).
%
% Initializes a halo orbit via its Jacobi constant CJAC. Alternative to the
% routine INIT_ORBIT(CR3BP, LI, TYPE, AFDIM, FAMILY, CST), which can
% initialize a halo orbit by an estimate of its vertical extension Az, via
% the parameter AFDIM.
%
% Inputs:
% 
% 1. CR3BP the structure containing the parent CR3BP
% 2. LI the structure containing the parent libration point.
% 3. FAMILY: the orbit family: either NORTHERN or SOUTHER for halo and
%    vertical lyapunov orbits. Always PLANAR for planar lyapunov orbits. 
%    Note that the family is forced to PLANAR for planar lyapunov orbits,
%    regardless of the input.
% 4. CJAC: the jacobi constant
% 5. CST the structure containing the numerical constants
%
% Outputs:
%
% 1. ORBIT the orbit structure , with the following fields:
%           - ORBIT.cr3bp       
%           - ORBIT.li          
%           - ORBIT.type        
%           - ORBIT.C           
%           - ORBIT.E           
%           - ORBIT.family      
%           - ORBIT.m           
%           - ORBIT.dm 
%           - ORBIT.status 
% See code for details on each of these fields.
%
% See also INIT_ORBIT
%
% BLB 2015

%CR3BP
orbit.cr3bp = cr3bp;

%Li
orbit.li = li;

%Family
orbit.family = family;

%Type
orbit.type = cst.orbit.type.HALO;

%Class
switch(family)
    case cst.orbit.family.NORTHERN
        orbit.m = 1;
    case cst.orbit.family.SOUTHERN
        orbit.m = 3;
end

%dm parameter (see Richardson)
switch li.number
    case 1
        orbit.dm = 2-orbit.m;
    case 2
        orbit.dm = orbit.m-2;  %BEWARE: Northern and Southern classes are 
                               %opposite to Class I and Class II 
                               %formulation (see Richardson 1980)
    case 3
        orbit.dm = 2-orbit.m;
    otherwise
        orbit.dm = NaN;  %NO dm if Li.number!=1,2,3
end

%Energy: 
% - orbit.C = Jacobi constant
% - orbit.E = -1/2*orbit.C = energy of the orbit
orbit.C = Cjac;
orbit.E = -0.5*Cjac;

%Status
orbit.status = cst.orbit.EMPTY;

end