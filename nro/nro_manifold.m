%--------------------------------------------------------------------------
% Plots an NRO orbit and its manifolds
%
% BLB 2016
%--------------------------------------------------------------------------
%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Change of parameters wrt default
default.plot.XY          = true;  %plot also the results in X-Y plane
default.plot.XZ          = true;  %plot also the results in X-Z plane
default.plot.YZ          = true;  %plot also the results in Y-Z plane

%% Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);

%% Orbit
%Initialize NRO
nro = init_nro(cr3bp, cr3bp.l2, cst.orbit.family.SOUTHERN, cst);

% Interpolation and plot
nro = nro_interpolation(cr3bp, nro, nro_init_EML2, default, cst, 'altitudeOfPerigee', 8000);


%% Manifold

% Event
event = init_event(cst.manifold.event.type.Y_SECTION,...                  %the event is triggered at a given angle...
                   0.0,...                                           %given in user data...
                   cst.manifold.event.isterminal.YES,...             %the trajectory stops at the first ocurrence...
                   cst.manifold.event.direction.ALL,...              %all direction are considered...
                   cr3bp.m2.pos, cst);                               %the center for the computation of the angle is the Moon
                    
% Manifold computation
% Unstable manifold
%-------------------------
manifold_branch_unstable_exterior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.EXTERIOR,...
                                          event);
                                      
manifold_branch_unstable_interior  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.INTERIOR,...
                                          event);

%Stable manifold
%-------------------------
manifold_branch_stable_exterior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.EXTERIOR,...
                                          event);
                                      
manifold_branch_stable_interior  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.INTERIOR,...
                                          event);
                                   
                                      
% Computation and plot
% Integration duration
t = 10;

%%
%Unstable interior
for theta = 0:0.1:1 %position on the orbit
    manifold_branch_unstable_interior = manifold_branch_computation(cr3bp, nro, manifold_branch_unstable_interior, theta, t, default, cst);
end

% % Unstable exterior
% for theta = 0:0.1:1 %position on the orbit
%     manifold_branch_unstable_exterior = manifold_branch_computation(cr3bp, nro, manifold_branch_unstable_exterior, theta, t, default, cst);
% end
% 
% % Stable exterior
% for theta = 0:0.1:1 %position on the orbit
%     manifold_branch_stable_exterior = manifold_branch_computation(cr3bp, nro, manifold_branch_stable_exterior, theta, t, default, cst);
% end
% 
% % Stable interior
% for theta = 0:0.1:1 %position on the orbit
%     manifold_branch_stable_interior = manifold_branch_computation(cr3bp, nro, manifold_branch_stable_interior, theta, t, default, cst);
% end