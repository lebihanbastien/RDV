%--------------------------------------------------------------------------
% Continuation example n°3: Continuation procedure to produce a discrete set within 
% the family of vertical lyapunov orbits of EML2
%
% BLB 2015
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
default.plot.XY              = false;     %plot also the results in X-Z plane

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit_1 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 16000, cst);
orbit_2 = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.VLYAP, cst.orbit.family.NORTHERN, 20000, cst);

%% Orbit computation
orbit_1 = orbit_computation(cr3bp, orbit_1, default, cst);
orbit_2 = orbit_computation(cr3bp, orbit_2, default, cst);

%% Continuation
%Waitbar
h = waitbar(0,'Computation in progress...');
%Loop
maxIter = 30;
for i = 1:maxIter
    yv = orbit_2.y0 + (orbit_2.y0 - orbit_1.y0);
    orbit_1 = orbit_2;
    orbit_2 = orbit_refinement(cr3bp, orbit_2, default, yv, cst);
    waitbar(i / maxIter);
end
close(h)

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end
