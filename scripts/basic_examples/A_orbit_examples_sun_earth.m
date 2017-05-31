%--------------------------------------------------------------------------
% This matlab file allows to build:
% - an EML2 planar lyapunov orbit (~90 000 km of radius)
% - an EML2 halo orbit (~90 000km of vertical extension)
% - an EML2 vertical lyapunov orbit (~90 000km of vertical extension)
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters

% 1. Decomment the next line to use only MATLAB routines (very slow!)
%--------------------------------------------------------------------------
%default.computation.type = cst.computation.MATLAB;

% 2. Decomment the next lines to change the absolute and relative precision during integration with MATLAB routines  (ode45)
%--------------------------------------------------------------------------
% default.ode45.AbsTol = 1e-12;
% default.ode45.RelTol = 1e-12;

% 3. Decomment the next lines to change the plotting options
%--------------------------------------------------------------------------
%default.plot.XZ        = true; %plot also the results in X-Z plane
%default.plot.YZ        = true; %plot also the results in Y-Z plane
%default.plot.diff_corr = true; %plot the differential correction steps

% 4. See parameters_default_init.m to see other options
%--------------------------------------------------------------------------

%% Environment init
cr3bp = init_CR3BP('SUN', 'EARTH', default);

%% Orbit init & computation for a planar lyapunov orbit
%Initialization
plyap = init_orbit(cr3bp, ...      % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.PLYAP, ...      % Planar lyapunov orbit
    cst.orbit.family.PLANAR, ...   % Planar class (useless here, since it is a planar lyapunov orbit
    90000, ...                      % Of Ax extension ~ 8000 km
    cst);                          % Numerical constants

%Computation
plyap = orbit_computation(cr3bp, plyap, default, cst);

%% Same for a halo orbit
%Initialization
halo = init_orbit(cr3bp, ...       % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.HALO, ...       % HALO orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    90000, ...                     % Of vertical extension ~ 10000 km
    cst);                          % Numerical constants

%Computation
halo = orbit_computation(cr3bp, halo, default, cst);

%% Same for a vertical lyapunov orbit
%Initialization
vlyap = init_orbit(cr3bp, ...      % Parent CR3BP
    cr3bp.l2, ...                  % Parent libration point
    cst.orbit.type.VLYAP, ...      % Vertical lyapunov orbit
    cst.orbit.family.NORTHERN, ... % Northern class
    90000, ...                     % Of vertical extension ~ 30000 km
    cst);                          % Numerical constants
%Computation
vlyap = orbit_computation(cr3bp, vlyap, default, cst);

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
    figure(4);
    view([-47 28]);
end
