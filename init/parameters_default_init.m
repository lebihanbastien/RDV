function default = parameters_default_init(cst)
% PARAMETERS_DEFAULT defines various default values for parameters used
% throughout the computations.
%
% DEFAULT = PARAMETERS_DEFAULT(CST) tnitializes various useful parameters 
% with default values, using constants nested in the structure CST.
%       
% BLB 2016

%--------------------------------------------------------------------------
%Integration
%--------------------------------------------------------------------------
default.ode45.RelTol  = 3e-14; %minimum allowed by ode45
default.ode45.AbsTol  = 1e-14; 

default.ode113.RelTol = 3e-14; %minimum allowed by ode113
default.ode113.AbsTol = 1e-14; 

default.ode87.RelTol  = 3e-14; 
default.ode87.AbsTol  = 1e-14; 

%--------------------------------------------------------------------------
%Differential correction
%--------------------------------------------------------------------------
default.diff_corr.precision = 1e-14;
default.diff_corr.type = cst.corr.Z0_FIXED;  %Rq: Z0_FIXED should be used with small Az, X0_FIXED otherwise. See orbit_refinement.m for details.
%[Deprecated]
default.diff_corr.isON = false;

%--------------------------------------------------------------------------
%Precision in libration points computation
%--------------------------------------------------------------------------
default.libp.precision = 1e-12;

%--------------------------------------------------------------------------
%Plotting or not?
%--------------------------------------------------------------------------
default.plot.orbit           = true;  %Global switch for plotting
default.plot.diff_corr       = false; %during differential correction
default.plot.XY              = true;  %Plot the XY view
default.plot.YZ              = false; %Plot the YZ view
default.plot.XZ              = false; %Plot the XZ view
default.plot.TD              = true;  %Plot the 3D view
default.plot.manifold_branch = true;  %during manifold computation
default.plot.LineSmoothing   = 'off';     %during manifold computation

default.plot.firstPrimDisp   = false;  %is the first primary (e.g. the Sun in the Sun-Earth system) displayed?
default.plot.allLibPoints    = false;  %are all libration points displayed?
default.plot.names           = false;  %are the names displayed?
default.plot.tdAxes          = false;  %are the pretty 3D axes displayed?
default.plot.bigPrimFac      = 1.0;        %the primaries appear bigPrimFac x bigger than they actually are (easier to see on screen)

%[Deprecated]
default.plot.halo_orbit      = true;  %during halo orbit computation (deprecated)
default.plot.lyap_orbit      = true;  %during lyap orbit computation (deprecated)

%--------------------------------------------------------------------------
% Computation type:
%   - cst.computation.MATLAB: only MATLAB routines
%   - cst.computation.MEX: MEX (compiled C) routines are included for
%   speed. Some functionalities are lost, such as user-defined event
%   routines.
%--------------------------------------------------------------------------
default.computation.type = cst.computation.MATLAB; %use MATLAB routines by default

end

