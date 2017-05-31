%--------------------------------------------------------------------------
% This matlab file allows to build an EML2 planar lyapunov orbit.
% Then, its exterior stable and unstable manifolds are computed
% and integrated up to given termination conditions:
% - The stable manifold is stopped at a given angle with respect to the
% Earth-Moon line, centered at the Earth.
% - The unstable manifold is stopped at x = x0.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
default.plot.TD = false; %no 3D plot is necessary for planar orbits

%% User input data
user.Ax = 8000;     %size of the orbit in [km]

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Orbit
orbit = init_orbit(cr3bp, ...                     % Parent CR3BP 
                   cr3bp.l2, ...                  % Parent libration point 
                   cst.orbit.type.PLYAP, ...      % Planar lyapunov orbit  
                   cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit 
                   user.Ax, ...                   % Of Ax extension ~ user.Ax
                   cst);                          % Numerical constants

%% Orbit computation
orbit = orbit_computation(cr3bp, orbit, default, cst);

%% Events
stable.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...        %the event is triggered at a given angle...
                          degtorad(-90),...                                %given in user data...
                          cst.manifold.event.isterminal.YES,...            %the trajectory stops at the first ocurrence...
                          cst.manifold.event.direction.INCREASING,...      %increasing directions are considered...
                          cr3bp.m1.pos, cst);                              %the center for the computation of the angle is the Earth
               
unstable.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...      %the event is triggered at a given value x = x0
                            degtorad(-90),...                              %given in user data...
                            cst.manifold.event.isterminal.YES,...          %the trajectory stops at the first ocurrence...
                            cst.manifold.event.direction.ALL,...           %all directions are considered...
                            cr3bp.m1.pos, cst);                            %the center for the computation of the angle is the Earth               

%% Manifold initialization                             
% The stable manifold is stopped when stable.event occurs
manifold_branch_stable = init_manifold_branch_event(cst.manifold.STABLE, ...
                                                    cst.manifold.EXTERIOR,...
                                                    stable.event);
                                                
% The unstable manifold is stopped when unstable.event occurs
manifold_branch_unstable = init_manifold_branch_event(cst.manifold.UNSTABLE, ...
                                                      cst.manifold.EXTERIOR,...
                                                      unstable.event);

%% Manifold computation
t = 30;

%% Stable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, theta, t, default, cst);
end


%% Unstable
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable = manifold_branch_computation(cr3bp, orbit, manifold_branch_unstable, theta, t, default, cst);
end