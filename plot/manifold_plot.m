function [] = manifold_plot(yv, orbit, manifold_branch, params, cst, varargin)
% MANIFOLD_PLOT plots a manifold leg on predefined figures.
%
% MANIFOLD_PLOT(YV, ORBIT, MANIFOLD_BRANCH, PARAMS, CST) plots the manifold
% leg MANIFOLD_BRANCH associated to the orbit ORBIT, on the plots required
% by the user in PARAMS, via the data stored in YV. This routines makes use
% of numerical constants defined in CST. By default:
%    - a stable manifold is plot in green, an unstable one in red. 
%    - all objects are plotted in SI units (km).
%    - the figure numbers are predefined in the following order:
%           -- figure(1): XY plot.
%           -- figure(2): XZ plot.
%           -- figure(3): YZ plot.
%           -- figure(4): XYZ (3D) plot.
%
% MANIFOLD_PLOT(YV, ORBIT, MANIFOLD_BRANCH, PARAMS, CST, COLOR) does the
% same thing but with a user-defined color COLOR.
%
% See also INITPLOT2D, INITPLOT3D
%
% BLB 2016.

%--------------------------------------------------------------------------
% Switch on the number of inputs to define the color
%--------------------------------------------------------------------------
switch(nargin)
    case 5
        %------------------------------------------------------------------
        % Color is function of the stability
        %------------------------------------------------------------------
        if(manifold_branch.stability == cst.manifold.STABLE)
            if(manifold_branch.way == cst.manifold.EXTERIOR)
                color = rgb('super green');%[51 102 0]/255;
            else
                color = rgb('super green');
            end
        else
            if(manifold_branch.way == cst.manifold.EXTERIOR)
                color = rgb('light red');%[153 0 0]/255;
            else
                color = rgb('light red');
            end
        end
    case 6 %a specific color has been provided
        color = varargin{1};     
    otherwise
        error('Wrong number of inputs.');
end


%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
if(params.plot.XY)
    manifold_plot_view(yv, 1, 2, orbit, params, color, 1);
end

if(params.plot.XZ)
    manifold_plot_view(yv, 1, 3, orbit, params, color, 2);
end

if(params.plot.YZ)
    manifold_plot_view(yv, 2, 3, orbit, params, color, 3);
end

if(params.plot.TD)
    manifold_plot_3D(yv, orbit, params, color, 4);
end

end

%--------------------------------------------------------------------------
% Plot a manifold on a given plane (XY, YZ or XZ). Smallest primary
% included.
%--------------------------------------------------------------------------
function [] = manifold_plot_view(yv, ip, jp, orbit, params, color, index)
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

%----------
%Plot
%----------
figure(index);
hold on
%Orbit
plot(yv(:,ip)*Lf,yv(:,jp)*Lf, 'Color',  color, 'LineWidth', 1.5);
end

%--------------------------------------------------------------------------
% Plot a manifold in 3D. Primaries included.
%--------------------------------------------------------------------------
function [] = manifold_plot_3D(yv, orbit, params,  color, index)
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

%----------
%Manifold branch
%----------
figure(index);
hold on
grid on
%Manifold branch
plot3(yv(:,1)*Lf, yv(:,2)*Lf,yv(:,3)*Lf, 'Color', color, 'LineWidth', 1.5);
    
end