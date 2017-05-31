function [Xtilde,t] = orbit_postprocess_modified_loic(cr3bp, orbit_target,orbit_chaser, y0target, y0chaser, params, cst)
%   ORBIT_POSTPROCESS Postprocessing routine for generic orbit.
%
%   ORBIT_POSTPROCESS(CR3BP, ORBIT, PARAMS, CST) computes some key elements
%   of the orbit. Namely:
%
%   - Either the couple (Az, Azdim) - vertical extension for halo and
%   vertical orbits, or the couple (Ax, Axdim), maximum planar extension
%   for planar lyapunov orbits.
%   - orbit.C: the jacobian constant
%   - orbit.E: the energy
%   - orbit.yv: the state along the orbit on a given grid over the interval
%   [0, orbit.T].
%   - orbit.monodromy: the monodromy matrix.
%   - orbit.eigenvalues: the eigenvalues of the monodromy matrix in vector
%   form.
%   - orbit.stable_direction: the stable eigenvector of the monodromy
%   matrix.
%   - orbit.unstable_direction: the unstable eigenvector of the monodromy
%   matrix.
%
% See Koon et al. 2006, chapter 6, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
%   This routine makes the following field of the orbit structure:
%   - orbit.y0:         initial conditions.
%   - orbit.T12:        half period.
%   - orbit.T:          period.
%
%   BLB 2016

%--------------------------------------------------------------------------
% True Az/Ax
%--------------------------------------------------------------------------
%Integration over one half orbit
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
    [~,y_target] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu), [0 orbit_target.T12], orbit_target.y0(1:6), options);
    yve = y_target(end,:);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    [~, yve] = ode78_cr3bp(0.0, orbit_target.T12, orbit_target.y0(1:6), 6, cr3bp.mu);
end

% Choose the maximum value between yve and yv0.
switch(orbit_target.type)
    case {cst.orbit.type.HALO, cst.orbit.type.VLYAP}
        orbit_target.Az    = max(abs(yve(3)), abs(orbit_target.y0(3)));
        orbit_target.Azdim = orbit_target.Az*cr3bp.L;
    case cst.orbit.type.PLYAP
        orbit_target.Ax    = max(abs(yve(1)), abs(orbit_target.y0(1)));
        orbit_target.Axdim = orbit_target.Ax*cr3bp.L;
end

%--------------------------------------------------------------------------
% Energy
%--------------------------------------------------------------------------
orbit_target.C = jacobi(orbit_target.y0, cr3bp.mu);  %jacobi constant
orbit_target.E = -0.5*orbit_target.C;                %energy

%--------------------------------------------------------------------------
% Change the period for Vertical orbits
%--------------------------------------------------------------------------
if(strcmp(orbit_target.type, cst.orbit.type.VLYAP))
    orbit_target.T12 = 2*orbit_target.T12;
    orbit_target.T   = 2*orbit_target.T;
end

%--------------------------------------------------------------------------
%Integration over one orbit (6 variables: state )
%--------------------------------------------------------------------------
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
    
    % simulation de deux orbites nro 
    %[t_target, y_target] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 orbit_target.T],orbit_target.y0(1:6), options);  
    %[t_chaser, y_chaser] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 orbit_chaser.T],orbit_chaser.y0(1:6), options);
    
    %Ttilde = zeros(length(y_target(:,1)),1);
    %Xtilde = zeros(6,length(y_target(:,1)));
%     x1t = y_target(:,1);
%     y1t = y_target(:,2);
%     z1t = y_target(:,3);
    %A0 = Atilde(cr3bp,x1t,y1t,z1t)
    tf = orbit_target.T;
    % position initiale cible
    y0(1:6) = y0target;
    % position initiale chasseur
    y0(7:12) = y0chaser;
    % position initiale rho = cible-target 
    y0(13:18) = (y0(7:12)-y0(1:6));
    [t,Xtilde] = ode113(@(t,y)cr3bp_relative_6(t,y,cr3bp),[0,10*tf], y0, options);

    yf = y_target(end,:);
    
    
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    [~, yf, y_target] = ode78_cr3bp(0.0, orbit_target.T, orbit_target.y0, 42, cr3bp.mu);
    %[~, yf, tv, yv] = ode78_cr3bp(0.0, orbit.T, orbit.y0, 42, cr3bp.mu);
end

%--------------------------------------------------------------------------
% Save the state on a grid over one period
%--------------------------------------------------------------------------
orbit_target.yv = y_target;
%orbit_target.tv = t_target;

%--------------------------------------------------------------------------
% Max/Min distance with each primary
%--------------------------------------------------------------------------
[orbit_target.minDistToM1, ~, orbit_target.maxDistToM1]  = distToPrimary(y_target, cr3bp.m1);
[orbit_target.minDistToM2, minDistIndex, orbit_target.maxDistToM2, maxDistIndex] = distToPrimary(y_target, cr3bp.m2);

%--------------------------------------------------------------------------
% 'Perigee' wrt M2 (Moon in the Earth-Moon system)
%--------------------------------------------------------------------------
orbit_target.perigee.radius   = orbit_target.minDistToM2;
orbit_target.perigee.altitude = orbit_target.minDistToM2-cr3bp.m2.Rm/cr3bp.L;
orbit_target.perigee.position = y_target(minDistIndex, 1:3);

%--------------------------------------------------------------------------
% 'Apogee' wrt M2 (Moon in the Earth-Moon system)
%--------------------------------------------------------------------------
orbit_target.apogee.altitude = orbit_target.maxDistToM2;
orbit_target.apogee.position = y_target(maxDistIndex, 1:3);

%--------------------------------------------------------------------------
% Linear algebra (monodromy matrix, etc)
%--------------------------------------------------------------------------
%Monodromy matrix
% orbit_target.monodromy = eye(6);
% for i = 1 : 6
%     for j = 1 : 6
%         m = 6*(i-1) + j;
%         orbit_target.monodromy(i,j) = yf(m+6);
%     end
% end
% 
% %Eigen
% [V,E] = eig(orbit_target.monodromy);
% 
% %Eigenvalues
% for i = 1:6
%     orbit_target.eigenvalues(i) = E(i,i);
% end
% 
% %Stable and unstable direction (linear approx of the manifolds)
% [~, posEigen] = min(abs(orbit_target.eigenvalues));
% orbit_target.stable_direction = V(:,posEigen);
% [~, posEigen] = max(abs(orbit_target.eigenvalues));
% orbit_target.unstable_direction = V(:,posEigen);
% 
end


function [minDist, minDistIndex, maxDist, maxDistIndex] = distToPrimary(yv, primary)
% [MINDIST, MINDISTINDEX, MAXDIST, MAXDISTINDEX] = DISTTOPRIMARY(YV, PRIMARY)
% computes the the min/max distance between the current state
% YV and the center of the primary PRIMARY.
%
% BLB 2016

% Get the difference in position for each values in yv
ydist  = bsxfun(@minus,yv(:,1:3), primary.pos);
% Get the norm of each difference vector
yrdist = arrayfun(@(idx) norm(ydist(idx,:)), 1:size(ydist,1));
% Get the mininum
[minDist, minDistIndex] = min(yrdist);
% Get the maximum
[maxDist, maxDistIndex] = max(yrdist);
end