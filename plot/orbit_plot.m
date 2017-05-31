function [] = orbit_plot(orbit, params, varargin)
% ORBIT_PLOT plots an orbit on predefined figures.
%
% ORBIT_PLOT(ORBIT, PARAMS) plots the orbit ORBIT, on the plots required
% by the user in PARAMS. By default:
%    - an orbit is plot in blue. 
%    - all objects are plotted in SI units (km).
%    - the figure numbers are predefined in the following order:
%           -- figure(1): XY plot.
%           -- figure(2): XZ plot.
%           -- figure(3): YZ plot.
%           -- figure(4): XYZ (3D) plot.
%
% ORBIT_PLOT(ORBIT, PARAMS, COLOR) does the same thing but with a 
% user-defined color COLOR.
%
% See also INITPLOT2D, INITPLOT3D
%
% BLB 2016.

%--------------------------------------------------------------------------
% Switch on the number of inputs
%--------------------------------------------------------------------------
switch(nargin)
    case 2
        color = rgb('dark blue');    
    case 3 %a specific color has been provided
        color = varargin{1};     
    otherwise
        error('Wrong number of inputs.');
end

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------
if(params.plot.XY)
    halo_orbit_plot_view(orbit, 1, 2, 1, params, color);
end

if(params.plot.XZ)
    halo_orbit_plot_view(orbit, 1, 3, 2, params, color);
end

if(params.plot.YZ)
    halo_orbit_plot_view(orbit, 2, 3, 3, params, color);
end

if(params.plot.TD)
    halo_orbit_plot_3D(orbit, 4, params, color);
end

end

%--------------------------------------------------------------------------
% Plot an orbit on a given plane (XY, YZ or XZ). Smallest primary
% included
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_view(orbit, ip, jp, index, params, varargin)
%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%Factor
%----------
Lf = cr3bp.L;

%----------
%Plot
%----------
%If the figure did not exist before, we set the basic environment.
if(~ishandle(index))
    initplot2D(index, cr3bp, params, orbit.li, ip, jp);
end

% ----------
% Orbit
% ----------
figure(index);
hold on;
switch(nargin)
    case 5
        plot(orbit.yv(:,ip)*Lf,orbit.yv(:,jp)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);
    case 6
        plot(orbit.yv(:,ip)*Lf,orbit.yv(:,jp)*Lf, 'Color',  varargin{1}, 'LineWidth', 1.5);
    otherwise
        error('Wrong number of inputs.');
end

end

%--------------------------------------------------------------------------
% Plot an orbit in 3D. Various options are available through params
%--------------------------------------------------------------------------
function [] = halo_orbit_plot_3D(orbit, index, params, varargin)

%----------
%Cr3bp
%----------
cr3bp = orbit.cr3bp;

%----------
%Factor
%----------
Lf = cr3bp.L;

%----------
%Plot
%----------
%If the figure did not exist before, we set the basic environment.
if(~ishandle(index))
    initplot3D(index, cr3bp, params);
end

% ----------
% Orbit
% ----------
figure(index);
hold on;
switch(nargin)
    case 3
        plot3(orbit.yv(:,1)*Lf, orbit.yv(:,2)*Lf,orbit.yv(:,3)*Lf, 'Color',  rgb('dark blue'), 'LineWidth', 1.5);
    case 4
        plot3(orbit.yv(:,1)*Lf, orbit.yv(:,2)*Lf,orbit.yv(:,3)*Lf, 'Color',  varargin{1}, 'LineWidth', 1.5);
    otherwise
        error('Wrong number of inputs.');
end


figure(index);
hold on;


end

