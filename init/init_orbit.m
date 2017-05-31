function orbit =  init_orbit(cr3bp, li, type, family, Afdim, cst)
% INIT_ORBIT initializes an orbit (halo, vertical of planar lyapunov)
% of the CRTBP.
%
% Inputs:
% 
% - cr3bp the structure containing the parent CR3BP
% - li the structure containing the parent libration point.
% - type the orbit type: either Halo, vertical lyapunov, planar lyapunov
% - family: the orbit family: either NORTHERN or SOUTHER for halo and
%    vertical lyapunov orbits. Always PLANAR for planar lyapunov orbits. 
%    Note that the family is forced to PLANAR for planar lyapunov orbits,
%    regardless of the input.
% - Afdim: the size of the orbit (in km). Equal to the vertical extension
%    Azdim in the case of halo and vertical lyapunov orbits; equal to the
%    planar extension Axdim in the case of a planar lyapunov orbit.
% - cst the structure containing the numerical constants
%
% Outputs:
%
% - orbit: the structure orbit, with the following fields:
%           - orbit.cr3bp       |
%           - orbit.li          |
%           - orbit.type        |
%           - orbit.Azdim       |   for Halo/Vertical lyapunov
%           - orbit.Az          |
%           - orbit.family      |
%           - orbit.m           |
%           - orbit.dm          |
%           - orbit.status      |
%    and  
%           - orbit.cr3bp       |
%           - orbit.li          |
%           - orbit.type        |
%           - orbit.Axdim       |   for Planar lyapunov
%           - orbit.Ax          |
%           - orbit.family      |
%           - orbit.status      |   
%
% BLB 2016


%CR3BP
orbit.cr3bp = cr3bp;

%Li
orbit.li = li;

%Type
orbit.type  = type;

%Family
switch(type)
    case {cst.orbit.type.HALO, cst.orbit.type.VLYAP}
        orbit.family = family;
    case cst.orbit.type.PLYAP
        %family forced to PLANAR, regardless of the input
        orbit.family = cst.orbit.family.PLANAR; 
end


%Class
switch(family)
    case cst.orbit.family.NORTHERN
        orbit.m = 1;
    case cst.orbit.family.SOUTHERN
        orbit.m = 3;
    otherwise
        %Class does not exist for planar orbits, NaN by default
        orbit.m = 0;
end

%dm parameters (see Richardson)
switch(li.number)
    case 1
        orbit.dm = 2-orbit.m;
    case 2
        %BEWARE: Northern and Southern classes
        %are opposite to Class I and Class II
        %formulation (see Richardson 1980)
        orbit.dm = orbit.m-2;
    case 3
        orbit.dm = 2-orbit.m;
    otherwise
        orbit.dm = 0;  %NO dm if Li.number!=1,2,3
end


%Extension
switch(type)
    case {cst.orbit.type.HALO, cst.orbit.type.VLYAP}
        %--------------
        %Afdim = Azdim
        %--------------
        orbit.Azdim = Afdim;
        orbit.Az = Afdim/cr3bp.L;
        %Az estimate
        orbit.Azdim_estimate = Afdim;
        orbit.Az_estimate = Afdim/cr3bp.L;

    case cst.orbit.type.PLYAP
        %--------------
        %Afdim = Axdim
        %--------------
        %Ax
        orbit.Axdim = Afdim;
        orbit.Ax = Afdim/cr3bp.L;
        %Ax estimate
        orbit.Axdim_estimate = Afdim;
        orbit.Ax_estimate = Afdim/cr3bp.L;
end

%Status
orbit.status = cst.orbit.EMPTY;

end