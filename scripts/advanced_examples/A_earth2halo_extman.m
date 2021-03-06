%--------------------------------------------------------------------------
% Advanced example n°1: 
%  
% This matlab file compute an Earth-to-Halo (EML2) transfer with:
% - A first maneuver @LEO
% - A second maneuver at the injection point in an exterior stable manifold
%
% The transfer is computed for several points
% on the halo orbit and, in the corresponding loop, the previous DeltaV 
% is used to compute the transfer for the next point, 
% providing that the previous solution is valid
%
% Both maneuvers are computed backwards in time, starting from a given
% point on an EML2-halo orbit. The LEO is defined only by its altitude
% (the final output.earth.inclination is governed by the correction scheme)
%
% The algorithm of the differential correction prodedure @LEO has been
% inpired by Gordon 2008 (Master thesis)
%
% This matlab file makes of the abacus
% ./data/halo_init_matrix_EML2.dat to generate an EML2 halo orbit and its
% exterior stable and unstable manifolds.
%
% WARNING: this computation includes a first guess fo the deltaV of the 
% second maneuver (at the insertion point of the stable/unstable manifold 
% branch).
% This first guess is implemented from heuristic considerations, and has a
% great influence on the final result.
%
% BLB 2016
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
default.plot.TD = false;  %do not plot the 3D results.

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, 8000, cst);
%Interpolation matrix
halo_init = halo_init_EML2;

%% User input data
%--------------------------------------------------------------------------
%May be modified
%--------------------------------------------------------------------------
user.hLEO       = 185;                %Desired LEO altitude [km]
user.theta      = 0;                  %Arbitrary position on the orbit, in [0 1]
user.showSteps  = false;              %To show/hide the steps of the corrective scheme
user.maxGapV    = 5.0;                %Maximum velocity gap allowed at injection point [km/s]
user.tangentManeuver = false;         %Is the Flyby Maneuver forced to be tangent? Very restrictive! May lead to no solution
user.fbangle  = degtorad(45);         %Termination angle of the manifold wrt the Earth-Moon line

%--------------------------------------------------------------------------
%Computed from user inputs, or arbitrary
%--------------------------------------------------------------------------
user.hLEOa = user.hLEO/cr3bp.L;  %Desired LEO altitude [adim]
user.t0    = 10;                 %Integration duration (arbitrary, will not be reached because of the termination @lunar flyby
user.isBCP = false;              %Are we using the Bicircular Problem? No

%% Orbit computation
orbit = halo_orbit_interpolation(cr3bp, orbit, halo_init, default, cst);

%% Events
%At the Moon
moon.event = init_event(cst.manifold.event.type.ANGLE_SECTION,...          %the event is triggered at a given angle...
                        user.fbangle,...                                   %given in user data...
                        cst.manifold.event.isterminal.YES,...              %the trajectory stops at the first ocurrence...
                        cst.manifold.event.direction.ALL,...               %all direction are considered...
                        cr3bp.m1.pos, cst);                                %the center for the computation of the angle is the Earth
               
earth.event = init_event(cst.manifold.event.type.FLIGHT_PATH_ANGLE,...     %the event is triggered when the flight path angle is...
                         0.0,...                                           %equal to zero...
                         cst.manifold.event.isterminal.YES,...             %the trajectory stops at the first ocurrence...
                         cst.manifold.event.direction.INCREASING,...       %all direction are considered...
                         cr3bp.m1.pos, cst);                               %the center for the computation of the angle is the Earth

%% Additional options
%Associated ode options
earth.options = odeset('Events',@(t,y)odezero_flightpathangle(t,y,earth.event),...
                       'Reltol', default.ode45.RelTol,...
                       'Abstol', default.ode45.AbsTol);

% Position of the primaries
moon.position  = cr3bp.m2.pos;
earth.position = cr3bp.m1.pos;


%% Manifold initialization
manifold_branch_stable  = init_manifold_branch_event(cst.manifold.STABLE,...    %the manifold is taken stable...
                                                     cst.manifold.EXTERIOR,...  %and interior...
                                                     moon.event);               %the integration is stopped when a certain angle wrt to the moon is reached

%% Loop on the position in L2
%Waitbar
h = waitbar(0,'Computation in progress...');
%--------------------------------------------------------------------------
% A rough first guess is used in lfb for the first entry in the loop.
% Then, output.manifold.deltaV can be used as first guess for the next
% iteration.
%--------------------------------------------------------------------------
isPreviousSolution = false;

% For saving outputs throughout the loop
it = 1;
for theta = 0:0.1:1
    %----------------------------------------------------------------------
    %New starting point
    %----------------------------------------------------------------------
    user.theta = theta;
    
    %----------------------------------------------------------------------
    %Manifold computation
    %----------------------------------------------------------------------
    manifold_branch_stable = manifold_branch_computation(cr3bp, orbit, manifold_branch_stable, user.theta, user.t0, default, cst);
    
    %----------------------------------------------------------------------
    % LFB
    %----------------------------------------------------------------------
    if(isPreviousSolution) % a first guess is added from a previous solution
        [output, isPreviousSolution] = blt(manifold_branch_stable, cr3bp, earth, moon, user, default, cst, output.manifold.deltaV);
    else % no previous solution
        [output, isPreviousSolution] = blt(manifold_branch_stable, cr3bp, earth, moon, user, default, cst); 
    end
    
    %----------------------------------------------------------------------
    % Saved outputs throughout the loop
    %----------------------------------------------------------------------
    if(isPreviousSolution)
        deltaVD(it) = output.deltaV_dim;
        thetaV(it)  = user.theta;
        Ttot(it)    = output.Ttot_dim; 
        it = it+1;
    end
    
    %----------------------------------------------------------------------
    %Uncomment this line if you do not want to use the first guess from a
    %previous solution
    %----------------------------------------------------------------------
    %isPreviousSolution = false;
    
    %----------------------------------------------------------------------
    % Waitbar
    %----------------------------------------------------------------------
    waitbar(theta);
end
close(h)


%% Maneuver cost vs position on the Halo orbit
figure;
hold on
grid on
xlabel('Position on the Halo orbit in [0 1]');
ylabel('Total maneuver cost [km/s]');
plot(thetaV, deltaVD, 'ob', 'MarkerFaceColor', 'black'); 


%% Maneuver cost vs time of flight
figure;
hold on
grid on
xlabel('Time of flight [days]');
ylabel('Total maneuver cost [km/s]');
plot(Ttot, deltaVD, 'ob', 'MarkerFaceColor', 'black'); 

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end
