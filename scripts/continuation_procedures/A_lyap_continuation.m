%--------------------------------------------------------------------------
% Continuation example n°1: Continuation procedure to produce a discrete 
% set within the family of planar lyapunov orbits of EML2
%
% BLB 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 1000, cst);
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.PLYAP, cst.orbit.family.PLANAR, 3000, cst);

%% Orbit computation
orbit_1 = orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation
maxIter = 23;
%Waitbar
h = waitbar(0,'Computation in progress...');
%Loop
for i = 1:maxIter
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    waitbar(i / maxIter);
end
% Close waitbar
close(h)

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end

