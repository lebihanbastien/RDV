%--------------------------------------------------------------------------
% Example n°6: Comparison of two different ways to compute halo orbits:
%
%   - The interpolation from a given abacus (halo_orbit_interpolation).
%   - The computation from a third-order first guess (halo_orbit_computation).
%
% This matlab file does two things:
%
%  1. It makes us of the abacus
% ./data/halo_init_matrix_EML1.dat and ./data/halo_init_matrix_EML2.dat
% to generate two halo orbits around L1 and L2 of same Az (15000km)
% via the routine:
%                   halo_orbit_interpolation
%  2. It computes the orbit of same Az_estimate starting from the 3rd
% order approximation of Richardson, via the routine:
%                   halo_orbit_computation
%--------------------------------------------------------------------------
% Important remarks:
%
% A. If MATLAB routines are used:
% When the precision of the integration is high (see default.ode45
% structure), the difference in CPU time between the interpolation process 
% and the computation from an initial 3rd order guess can be important.
%
% B. For the Earth-Moon system, one can try to set Az = 35000 km. 
% For this value, the third order approximation of Richardson is clearly
% not valid anymore: the difference between the interpolated orbits (the
% good ones) and the orbits computed from the 3rd order approx is huge.
% For even higher values of Az, the computation based on Richardson's 3rd
% order diverges.
%
% C. However, the interpolation process requires abacuses such as 
% ./data/halo_init_matrix_EML1.dat in order to generate the solutions.
% in these version of the code, only the Earth-Moon L1,2 abacuses are
% provided. For other computation (either different points or different
% systems), the halo_orbit_computation routine MUST be used, keeping in
% mind that it is valid for only a certain range of Az values.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------

%% Initialization: reboot, addpath, constants, default parameters. See init.m
init;

%% Inner changes from default parameters
default.plot.orbit = false; % Plotting is at the end so that the CPU time 
                            % used for the plots is not taken into account 
                            % in the above computations

%% User data
user.Az = 30000; %Az in [km]

%% Structures init
%Environment
cr3bp = init_CR3BP('EARTH', 'MOON', default);
             
%% Orbit computation with third order approximation

%Init
hl2_comp = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, user.Az, cst);
hl1_comp = init_orbit(cr3bp, cr3bp.l1,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, user.Az, cst);    

tic;
hl2_comp = orbit_computation(cr3bp, hl2_comp, default, cst);
hl1_comp = orbit_computation(cr3bp, hl1_comp, default, cst);
T = toc;

% Conclusion
fprintf('End of computation via ORBIT_COMPUTATION in %5.5f s\n', T);
fprintf('EML1: Desired Az was %5.2fkm, actual is %5.2fkm \n', hl1_comp.Azdim_estimate, hl1_comp.Azdim);
fprintf('EML2: Desired Az was %5.2fkm, actual is %5.2fkm \n', hl2_comp.Azdim_estimate, hl2_comp.Azdim);
%% Orbit computation with interpolation

%Init
hl2_int = init_orbit(cr3bp, cr3bp.l2,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, user.Az, cst);
hl1_int = init_orbit(cr3bp, cr3bp.l1,  cst.orbit.type.HALO, cst.orbit.family.NORTHERN, user.Az, cst); 

tic;
hl2_int = halo_orbit_interpolation(cr3bp, hl2_int, halo_init_EML2, default, cst);
hl1_int = halo_orbit_interpolation(cr3bp, hl1_int, halo_init_EML1, default, cst);
T = toc;

% Conclusion
fprintf('End of computation via HALO_ORBIT_INTERPOLATION in %5.5f s\n', T);
fprintf('EML1: Desired Az was %5.2fkm, actual is %5.2fkm \n', hl1_int.Azdim_estimate, hl1_int.Azdim);
fprintf('EML2: Desired Az was %5.2fkm, actual is %5.2fkm \n', hl2_int.Azdim_estimate, hl2_int.Azdim);

%% Plotting is at the end so that the CPU time used for the plots is not taken into account in the above computations
orbit_plot(hl2_comp, default, rgb('dark red'));
orbit_plot(hl1_comp, default, rgb('dark red'));
orbit_plot(hl2_int, default, rgb('dark green'));
orbit_plot(hl1_int, default, rgb('dark green'));

%% Change the orientation of the 3D plot, if it exists
if(any(findall(0,'Type','Figure')==4))
   figure(4);
   view([-47 28]);
end