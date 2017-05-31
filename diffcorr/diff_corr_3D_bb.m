function orbit = diff_corr_3D_bb(v0, cr3bp , orbit, xi0, xif, params, cst)
% DIFF_CORR_3D_BB Differential correction to compute 3D periodic orbits of
% the CRTBP, symmetric with respect to the xz-plane (halo, vertical lyapunov).
%
% This routines iteratively correct the initial state v0 via an iterative
% Newton method. This method is illustrated below with the example of the
% Halo orbits.
%
% For Halo orbits, at each step:
% The free variables are
%       X0 = [x0 vy0]^T 
% or
%       X0 = [z0 vy0]^T
%
% The constraint vector is FX = F(X) = [y vx vz]^T. The goal is to
% satisfy the constraint FX = 0, so as to produce a periodic orbit
% symmetric with respect to the xz-plane.
%
% Starting with initial condition Y0 = [x0 y0 z0 vx0 vy0 vz0]^T,
% the equations of motion of the CRTBP are integrated until the xz-plane is
% reached.
% With this choice the first constraint y = 0 is automatically satisfied.
%
% Then, the Newton's correction dX0 to be applied so as to obtain
%               FXr = [vx vz]^T = 0
% satisfies the following equation:
%
% Case (i):
%   FXr                  Af               ppf               Bf        dX0
% [dvx]  = ( [Phi41  Phi45] - 1/ypoint * [ax] * [Phi21   Phi25] ) * [dx0 ]
% [dvz]      [Phi61  Phi65]              [az]                       [dvy0]
%
% Case (ii):
%                      Af                 ppf           Bf
% [dvx]  = ( [Phi43  Phi45] - 1/ypoint * [ax] * [Phi23   Phi25] ) * [dz0 ]
% [dvz]      [Phi63  Phi65]              [az]                       [dvy0]
%
%  where Phi is the State Transition Matrix (STM), numerically
%  integrated along with the state.
%
% This correction is applied iteratively on the initial state until a
% use-defined precision threshold is reached.
% The vectors X0 and FXr as well as the equations (i) and (ii)
% are encoded in the inputs xi0 and xif. Example in case (i):
%           xi0 = [3, 5];  %X0  = [z0 vy0]^T is corrected
%           xif = [4, 6];  %FXr = [vx vz]^T = 0 is targeted
%
% The vectors X0 and FXr may differ from the ones written here for
% other types of orbits (e.g. vertical lyapunov orbits).
%
% See Howell 1984, for details <a href="matlab: 
% web('http://adsabs.harvard.edu/full/1984CeMec..32...53H','-browser')">(link)</a>.
%
%
% At the end of the routine, the following elements are updated:
% - orbit.y0:         initial conditions.
% - orbit.T12:        half period.
% - orbit.T:          period.
%
% BLB 2016.

%--------------------------------------------------------------------------
%Iteration counts
%--------------------------------------------------------------------------
iter = 0;

%--------------------------------------------------------------------------
%Type of event: reaching the y = 0 plane
%--------------------------------------------------------------------------
options = odeset('Events', @odezero_y,'Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
section = cst.manifold.event.type.Y_SECTION;
si = 2; % y = 0 is targeted.

%Associated event structure for MEX routines
val_par = init_event(section,...
                     0.0,...
                     cst.manifold.event.isterminal.YES,...
                     cst.manifold.event.direction.ALL,... 
                     cr3bp.m1.pos,...
                     cst);
                 
        
%--------------------------------------------------------------------------
% Differential correction loop
%--------------------------------------------------------------------------
while(true)
    iter = iter+1;
    
     if(iter > 50)
        disp('WARNING: maximum iterations reached in differential_correction');
        break;
     end
    
    %----------------------------------------------------------------------
    % Integration stops at y=0
    %----------------------------------------------------------------------
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        [~,yv,te,ve,~] = ode113(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 10],v0,options);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        if(params.plot.diff_corr)
            [te, ve, ~, yv] = ode78_cr3bp_event(0.0, 10, v0, 42, cr3bp.mu, val_par);
        else
            [te, ve] = ode78_cr3bp_event(0.0, 10, v0, 42, cr3bp.mu, val_par);
        end
        
    end
    
    
    %----------------------------------------------------------------------
    % Update the final state: FX = [vx vz]^T.
    %----------------------------------------------------------------------
    FXr = [ve(xif(1)) ; ve(xif(2))];
 
    %----------------------------------------------------------------------
    % Compute the first order correction with Newton's method
    %----------------------------------------------------------------------
    % Example, for HALO orbits, and corr_type == cst.corr.Z0_FIXED:
    %  FXr               Af                   ppf         Bf              dX0
    % [dvx]  = ( [Phi41  Phi45] - 1/ypoint * [ax] * [Phi21   Phi25] ) * [dx0 ]
    % [dvz]      [Phi61  Phi65]              [az]                       [dvy0]
    %
    %----------------------------------------------------------------------
    % Example, for HALO orbits, corr_type == cst.corr.X0_FIXED:
    %  FXr               Af                   ppf         Bf              dX0
    % [dvx]  = ( [Phi43  Phi45] - 1/ypoint * [ax] * [Phi23   Phi25] ) * [dz0 ]
    % [dvz]      [Phi63  Phi65]              [az]                       [dvy0]
    %
    %----------------------------------------------------------------------
    % Example, for VLYAP orbits, corr_type == cst.corr.Z0_FIXED:
    %  FXr               Af                   ppf         Bf              dX0
    % [dz]   = ( [Phi31  Phi35] - 1/ypoint * [zp]  * [Phi21   Phi25] ) * [dx0 ]
    % [dvx]      [Phi41  Phi45]              [ax]                        [dvy0]
    %
    %----------------------------------------------------------------------
    % Example, for VLYAP orbits, corr_type == cst.corr.X0_FIXED:
    %  FXr               Af                   ppf         Bf              X0
    % [dz]   = ( [Phi33  Phi35] - 1/ypoint * [zp]  * [Phi23   Phi25] ) * [dz0 ]
    % [dvx]      [Phi43  Phi45]              [ax]                        [dvy0]
    %----------------------------------------------------------------------
    %Af matrix
    Af = zeros(2);
    for i0 = 1:2
       for j0 = 1:2 
            Af(i0,j0) = ve(findIndix(xif(i0), xi0(j0)));
       end
    end
    % Bf matrix
    Bf = (1:2);
    Bf(1) = ve(findIndix(si,xi0(1)));
    Bf(2) = ve(findIndix(si,xi0(2)));
    
    %Derivative of ve at y=0 (t = te)
    vep = cr3bp_derivatives_6(te, ve, cr3bp.mu);
       
    %ppf vector
    ppf = (1:2)';
    ppf(1) = vep(xif(1))/ve(si+3);
    ppf(2) = vep(xif(2))/ve(si+3);
    
    %New Af
    Af = Af - ppf*Bf;

    %Stops if precision is good enough
    if(norm(FXr) < params.diff_corr.precision)
        break;
    end
    
    %----------------------------------------------------------------------
    % The first order correction is computed as dX0 = inv(Af)*FXr
    %----------------------------------------------------------------------
    dX0 = Af \ FXr;
    
    %----------------------------------------------------------------------
    %Updating initial state
    %----------------------------------------------------------------------
    v0(xi0(1)) = v0(xi0(1)) - dX0(1);
    v0(xi0(2)) = v0(xi0(2)) - dX0(2);   
    
    
    %----------------------------------------------------------------------
    % Plotting (potentially)
    %----------------------------------------------------------------------
    if(params.plot.diff_corr)
        plotDiffCorr(yv, cr3bp, iter);
    end
   
end

%Orbit update
orbit.y0  = v0;
orbit.T12 = te;   %1/2 period
orbit.T   = 2*te; %period

end

%--------------------------------------------------------------------------
% Return the right indix to find Phi(i0, j0) in ve
%--------------------------------------------------------------------------
function shiftedIndix = findIndix(i0, j0)
    shiftedIndix = 6 + 6*(i0-1)+j0;
end

%--------------------------------------------------------------------------
% Plotting the iterations
%--------------------------------------------------------------------------
function [] = plotDiffCorr(yv, cr3bp, iter)
    figure(1)
    hold on
    if(iter==1)
        plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'g');
    else
        plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'r');
    end
end