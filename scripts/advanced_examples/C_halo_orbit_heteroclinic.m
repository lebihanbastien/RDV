%--------------------------------------------------------------------------
% Advanced example n°3: 
%  
% This matlab file :
%
% 1. Makes of the abacus
% ./data/halo_init_matrix_EML2.dat to generate an EML2 halo orbit 
% with a given energy value
%
% 2. Makes of the abacus
% ./data/halo_init_matrix_EML1.dat to generate an EML1 halo orbit 
% with the same energy value
%
% 3. Starts a rough search for heteroclinic connections by computing the
% stable exterior manifold of the EML1 orbit and the unstable interior 
% of the EML2 orbit.
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
%Orbit
orbit_L1 = init_halo_orbit_energy(cr3bp, cr3bp.l1, cst.orbit.family.NORTHERN, 3.12, cst);
orbit_L2 = init_halo_orbit_energy(cr3bp, cr3bp.l2, cst.orbit.family.NORTHERN, 3.12, cst);


%% Orbit computation
orbit_L1 = halo_orbit_interpolation(cr3bp, orbit_L1, halo_init_EML1, default, cst);
orbit_L2 = halo_orbit_interpolation(cr3bp, orbit_L2, halo_init_EML2, default, cst);


%% Manifolds
manifold_L1  = init_manifold_branch_event(cst.manifold.STABLE,...
                                          cst.manifold.EXTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          1-cr3bp.mu,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
manifold_L2  = init_manifold_branch_event(cst.manifold.UNSTABLE,...
                                          cst.manifold.INTERIOR,...
                                          cst.manifold.event.type.X_SECTION,...
                                          1-cr3bp.mu,...
                                          cst.manifold.event.isterminal.YES,...
                                          cst.manifold.event.direction.ALL,...
                                          cr3bp.m1.pos,...
                                          cst);
                                      
%% Poincare sections
t = 10;
for theta = 0:0.05:1
    manifold_L1 = manifold_branch_computation(cr3bp, orbit_L1, manifold_L1, theta, t, default, cst);
    manifold_L2 = manifold_branch_computation(cr3bp, orbit_L2, manifold_L2, theta, t, default, cst);
end

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end