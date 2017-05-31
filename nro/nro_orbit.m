%--------------------------------------------------------------------------
% Plots an NRO orbit.
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = false;  %plot also the results in X-Y plane
default.plot.XZ          = false;  %plot also the results in X-Z plane
default.plot.YZ          = false;  %plot also the results in Y-Z plane
default.plot.TD          = true;   %plot also the results in three-D
%% Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit
%Initialize NRO
nro = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro = nro_interpolation(cr3bp, nro, nro_init_EML2, default, cst, 'altitudeOfPerigee', 500);

%% Plot T=f(z)
numberOfIteration = 80;
begin_altitude = 10;
end_altitude = 500;
Step = (end_altitude-begin_altitude)/numberOfIteration;
Time_period = zeros(numberOfIteration,2);

% Waitbar
h = waitbar(0,'Computation in progress...');

iteration = 1:numberOfIteration;

%Initialize NRO
nro_plot = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

for i = iteration
   altitude = begin_altitude + Step*(i-1);
   nro_plot = nro_interpolation(cr3bp, nro, nro_init_EML2, default, cst, 'altitudeOfPerigee', altitude);
   
   %Keep in memory
    % column 1 : Time period
    % cilumn 2 : Altitude
   Time_period(i,1) = nro_plot.T*cr3bp.T/(3600*24*2*pi);
   Time_period(i,2) = altitude;
   
   %----------------------------------------------------------------------
    % Waitbar
    %----------------------------------------------------------------------
    waitbar(i / length(iteration));
   
end

close(h)
close all

%Draw
plot(Time_period(:,2), Time_period(:,1));
xlabel('Altitude (km)');
ylabel('Period (day)');
title('Period = f(altitude)');
