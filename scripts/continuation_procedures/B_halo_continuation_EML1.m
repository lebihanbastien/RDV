%--------------------------------------------------------------------------
% This matlab script uses a pseudo-arclength continuation to built the
% Earth-Moon halo family about L1
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
% MEX files are used to fasten computation. Comment this line if MEX files
% do not work on you current version.
default.computation.type = cst.computation.MEX;

default.plot.XY          = true;  %plot also the results in X-Y plane
default.plot.XZ          = true;  %plot also the results in X-Z plane
default.plot.YZ          = true;  %plot also the results in Y-Z plane
default.plot.diff_corr   = false; %plot also the different steps in the differential correction procedures.
default.plot.orbit       = false;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%Initialize orbit
halo = init_orbit(cr3bp, ...        % Parent CR3BP
    cr3bp.l1, ...                   % Parent libration point
    cst.orbit.type.HALO, ...        % HALO orbit
    cst.orbit.family.NORTHERN, ...  % Northern class
    2000, ...                       % Of vertical extension ~ 10000 km
    cst);                           % Numerical constants

%% Pseudo-arclength continuation: initialization

%--------------------------------------------------------------------------
% First computation: a true orbit is computed with a simple 3-dimensionnal
% differential corrector scheme (no pseudo-arclength). The z component of
% the initial state is fixed.
%--------------------------------------------------------------------------
halo = orbit_computation(cr3bp, halo, default, cst, cst.corr.Z0_FIXED);

%--------------------------------------------------------------------------
% Second computation: initialized by the previous computation, a new orbit
% is obtained with a 4-dimensional differential corrector. The z component
% is now part of the free variables. The objects necessary to go on with
% the continuation are initialized in orbit.cont (free-variables vector,
% null vector of the Jacobian).
%--------------------------------------------------------------------------
halo = orbit_refinement(cr3bp, halo, default, halo.y0, cst, cst.corr.MIN_NORM);

%Arclength stepsize
halo.cont.ds = 0.005;

%% Pseudo-arclength continuation: loop

% Waitbar
h = waitbar(0,'Computation in progress...');

%--------------------------------------------------------------------------
% The pseudo-arclenght continuation is performed, on a given number of
% steps of size halo.cont.ds.
%--------------------------------------------------------------------------
% step vector
output.index = 1:400;
% Loop
for i = output.index
    %----------------------------------------------------------------------
    % Differential correction
    %----------------------------------------------------------------------
    halo = orbit_refinement_cont(cr3bp, halo, default, cst);
    
    %----------------------------------------------------------------------
    % Save outputs
    %----------------------------------------------------------------------
    output.orbit(i) = halo;
    
    %----------------------------------------------------------------------
    % Waitbar
    %----------------------------------------------------------------------
    waitbar(i / length(output.index));
end
close(h)

%% Postprocessing & plotting
% The following defining constraints for NRO are selected:
% - if the abscissa of the perigee with respect to the Moon projected on
%   the xy-plane ust intersect the Moon's radius.
% - The perigee altitude wrt to the Moon's surface must be positive.
%
% Select a given condition to isolate the NROs
for i = output.index
    %Get the perigee position
    pos = output.orbit(i).perigee.position;
    % If both conditions above are satisfied, we go on
    if(abs(pos(1) - cr3bp.m2.pos(1)) < cr3bp.m2.Rm/cr3bp.L && output.orbit(i).perigee.altitude > 0)
        output.isNRO(i) = true;
    else
        output.isNRO(i) = false;
    end
end

% Extract the different parameters in the vector output.orbit
output.C       = extractfield(output.orbit,'C');
output.T       = extractfield(output.orbit,'T');
output.Az      = extractfield(output.orbit,'Az');
output.Azdim      = extractfield(output.orbit,'Azdim');
output.perigee = cell2mat(extractfield(output.orbit,'perigee'));
output.altitudeOfPerigee = extractfield(output.perigee,'altitude');
output.perigeeDistanceToMoonCenter = extractfield(output.perigee,'radius');

%% Family of orbits
freq = 10;
nlength = length(output.index);

for i = 1:freq:nlength
    %Get the perigee position
    pos = output.orbit(i).perigee.position;
    % Orbit
    if(output.isNRO(i))
        orbit_plot(output.orbit(i), default, rgb('dark green'));
    else
        orbit_plot(output.orbit(i), default);
    end
    % Perigee position
    figure(4)
    hold on;
    plot3(pos(1)*cr3bp.L, pos(2)*cr3bp.L, pos(3)*cr3bp.L, 'ko', 'MarkerFaceColor', 'k');
end


%% Jacobi constant
figure;
hold on;
grid on;
plot(output.index, output.C, 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
plot(output.index(output.isNRO), output.C(output.isNRO), 'Color', rgb('dark green'), 'LineWidth', 2);
xlabel('Index [-]');
ylabel('Jacobi constant');
legend('halo', 'nro');

%% Altitude of the perigee
figure;
hold on;
grid on;
plot(output.index, output.altitudeOfPerigee*cr3bp.L, 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
plot(output.index(output.isNRO), output.altitudeOfPerigee(output.isNRO)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 2);
xlabel('Index [-]');
ylabel('Altitude of perigee [km]');
legend('halo', 'nro');

%% Altitude of the perigee vs Period
figure;
hold on;
grid on;
plot(output.T*cr3bp.T/(2*pi*86400), output.altitudeOfPerigee*cr3bp.L, 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
plot(output.T(output.isNRO)*cr3bp.T/(2*pi*86400), output.altitudeOfPerigee(output.isNRO)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 2);
xlabel('Period [days]');
ylabel('Altitude of perigee [km]');
legend('halo', 'nro');

%% Altitude of the perigee vs Az
figure;
hold on;
grid on;
plot(output.Azdim, output.altitudeOfPerigee*cr3bp.L, 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
plot(output.Azdim(output.isNRO), output.altitudeOfPerigee(output.isNRO)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 2);
xlabel('Az [km]');
ylabel('Altitude of perigee [km]');
legend('halo', 'nro');

%% Save the abacus
% Exctract the initial conditions
output.yv = vec2mat(extractfield(output.orbit,'y0'),42);

%--------------------------------------------------------------------------
% Build the abacus
%--------------------------------------------------------------------------
if(max(output.isNRO))
    % Initial conditions
    nro_init_EML1.initialConditions = output.yv(output.isNRO,1:6);
    
    % Altitude of perigee
    nro_init_EML1.altitudeOfPerigee = output.altitudeOfPerigee(output.isNRO);
    nro_init_EML1.altitudeOfPerigeeLimit(1) = min(nro_init_EML1.altitudeOfPerigee);
    nro_init_EML1.altitudeOfPerigeeLimit(2) = max(nro_init_EML1.altitudeOfPerigee);
    
    % Distance of the perigee wrt the center of the Moon
    nro_init_EML1.perigeeDistanceToMoonCenter = output.perigeeDistanceToMoonCenter(output.isNRO);
    nro_init_EML1.perigeeDistanceToMoonCenterLimit(1) = min(nro_init_EML1.perigeeDistanceToMoonCenter);
    nro_init_EML1.perigeeDistanceToMoonCenterLimit(2) = max(nro_init_EML1.perigeeDistanceToMoonCenter);
    
    % Vertical extension
    nro_init_EML1.Az = output.Az(output.isNRO);
    nro_init_EML1.AzLimit(1) = min(nro_init_EML1.Az);
    nro_init_EML1.AzLimit(2) = max(nro_init_EML1.Az);
    
    %--------------------------------------------------------------------------
    % Save abacus (only if necessary!)
    %--------------------------------------------------------------------------
    % save nro_init_EML1 nro_init_EML1
    
end

return;

%% Example of use

% Plot the orbit
default.plot.orbit = true;

%Initialize NRO
nro = init_nro(cr3bp, cr3bp.l1, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro = nro_interpolation(cr3bp, nro, nro_init_EML1, default, cst, 'Az', 80000);