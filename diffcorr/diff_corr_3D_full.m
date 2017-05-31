function orbit = diff_corr_3D_full(v0, cr3bp , orbit, params, cst)
% DIFF_CORR_3D_FULL Differential correction to compute for 3D periodic 
% orbits of the CRTBP symmetric, with respect to the xz-plane (halo, 
% vertical).
%
% This differential corrector is heavily based on the summary by Pavlak in
% his PhD thesis (2013), although the basic principles are much older than
% his work.
%
% See Pavlak 2013, section 3.3, for details <a href="matlab: 
% web('https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2013_Pavlak.pdf','-browser')">(link)</a>.
%
% This routines iteratively correct the initial state v0 via an iterative
% Newton method.
% Given the free variables X0:
%       X0 = [x0 z0 vy0 T12]^T, with vy0 = dot(y0), T12 = half period,
%
% given the targeted constraint FX = F(X) = [y vx vz]^T = 0,
%
% the correction dX0 to be applied to X0 satisfies the following equation:
%           FX = DFX*dX0 (1)
% with - FX = [y vx vz]^T,
%      - DFX the matrix defined as the Jacobian: DFX = dF(X)/dX0
%
% More precisely:
%
%   DFX = | Phi21   Phi23  Phi25  vy |
%	      | Phi41   Phi43  Phi45  ax |
%         | Phi61   Phi63  Phi65  az |
%
%   with vy = dot(y), ax = dot(vx) = dot(dot(x)), same for az.
%   and where Phi is the State Transition Matrix (STM), numerically
%   integrated along with the state.
%
%  The system (1), infradetermined. The minimum norm solution is chosen
%  here:
%           dX0 = DFX' * inv(DFX * DFX') * FX;
%
% At the end of the routine, the following elements are updated:
% - orbit.y0:         initial conditions.
% - orbit.T12:        half period.
% - orbit.T:          period.
% - orbit.cont.gv:    final vector of free variables for a potential continuation procedure.
% - orbit.cont.nv:    null vector of the Jacobian for a potential continuation procedure.
%
% BLB 2016.

%--------------------------------------------------------------------------
%Iteration counts
%--------------------------------------------------------------------------
iter = 0;

%--------------------------------------------------------------------------
% Get an estimate of the half period if it is not provided
%--------------------------------------------------------------------------
if(isfield(orbit, 'T12'))
    T12 = orbit.T12;
else
    %----------------------------------------------------------------------
    % If the half period is not provided from a previous computation, the
    % following code intent to compute an estimate of the half period by
    % integrating the equations of motion until the symmetry plane of the
    % orbit is reached.
    %
    % For example, for an halo orbit, the equations of motion should be
    % integrated until the plane y = 0 is reached.
    %----------------------------------------------------------------------
    disp( 'diff_corr_3D_full. No half period estimate was provided.');
    disp(['diff_corr_3D_full. An approximation of the half period ',...
        'is computed from a numerical integration until the ',...
        'symmetry plane of the orbit is reached.']);
    %----------------------------------------------------------------------
    % Numerical integration with ode113
    %----------------------------------------------------------------------
    options = odeset('Events', @odezero_y,'Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
    [~,~,T12] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),[0 10],v0(1:6),options);
end

%--------------------------------------------------------------------------
% Update the options: no event is necessary after this point.
%--------------------------------------------------------------------------
if(params.computation.type == cst.computation.MATLAB)
    options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
end

%--------------------------------------------------------------------------
% Differential correction loop
%--------------------------------------------------------------------------
while(true)
    iter = iter+1;
    %----------------------------------------------------------------------
    %Stops if too much iterations
    %----------------------------------------------------------------------
    if(iter > 50)
        disp('WARNING: maximum iterations reached in differential_correction');
        break;
    end
    
    %----------------------------------------------------------------------
    % Integration stops at t = T12
    %----------------------------------------------------------------------
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        [t,yv] = ode113(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 T12],v0,options);
        % Update final state
        te = t(end);
        ve = yv(end,:);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        if(params.plot.diff_corr)
            [te, ve, ~, yv] = ode78_cr3bp(0.0, T12, v0, 42, cr3bp.mu);
        else
            [te, ve] = ode78_cr3bp(0.0, T12, v0, 42, cr3bp.mu);
        end
    end
    
    %----------------------------------------------------------------------
    % Update the final state: FX = [y vx vz]^T.
    %----------------------------------------------------------------------
    FX = [ve(2) ; ve(4) ; ve(6)];
    
    %----------------------------------------------------------------------
    % Compute the first order correction with Newton's method:
    % Targeting the constraint FX = F(X) = 0, with the free variables X0,
    % The correction dX0 to applied to these variables satisfies the
    % equation:
    %           FX = DFX*dX0
    % With DFX the matrix defined as the Jacobian:
    %           DFX = dF(X)/dX0
    % In our case, the free variables are:
    %           X0 = [x0 z0 vy0 T12]^T
    %----------------------------------------------------------------------
    % Build the Jacobian DFX = dF(X)/dX0, with X0 = [x0, z0, vy0 T12]^T
    %
    %   DFX = | Phi21   Phi23  Phi25  vy |
    %	      | Phi41   Phi43  Phi45  ax |
    %         | Phi61   Phi63  Phi65  az |
    %
    %   with vy = dot(y), ax = dot(vx) = dot(dot(x)), same for az.
    %   and where Phi is the State Transition Matrix (STM).
    %----------------------------------------------------------------------
    DFX = zeros(3,4);
    
    % Concatenate the elements of the STM.
    DFX(1,1) = ve(findIndix(2,1));
    DFX(1,2) = ve(findIndix(2,3));
    DFX(1,3) = ve(findIndix(2,5));
    
    DFX(2,1) = ve(findIndix(4,1));
    DFX(2,2) = ve(findIndix(4,3));
    DFX(2,3) = ve(findIndix(4,5));
    
    DFX(3,1) = ve(findIndix(6,1));
    DFX(3,2) = ve(findIndix(6,3));
    DFX(3,3) = ve(findIndix(6,5));
    
    %Last column is [vy ax az]^T
    vep = cr3bp_derivatives_6(te, ve, cr3bp.mu);
    DFX(1,4) = vep(2);
    DFX(2,4) = vep(4);
    DFX(3,4) = vep(6);
    
    %----------------------------------------------------------------------
    % Stops if precision is good enough, but after the Jacobian has been
    % built, so that it is updated at the end.
    %----------------------------------------------------------------------
    if(norm(FX) < params.diff_corr.precision);
        break;
    end
    
    %----------------------------------------------------------------------
    % The first order correction is computed as the minimum norm solution
    %----------------------------------------------------------------------
    dX0 = DFX' * ( (DFX * DFX') \ FX);
    
    %----------------------------------------------------------------------
    % Updating the initial state
    %----------------------------------------------------------------------
    v0(1) = v0(1) - dX0(1);
    v0(3) = v0(3) - dX0(2);
    v0(5) = v0(5) - dX0(3);
    T12   = T12   - dX0(4);
    
    %----------------------------------------------------------------------
    % Plotting (potentially)
    %----------------------------------------------------------------------
    if(params.plot.diff_corr)
        plotDiffCorr(yv, cr3bp, iter);
    end
end

%--------------------------------------------------------------------------
% Orbit update
%--------------------------------------------------------------------------
orbit.y0  = v0;
orbit.T12 = T12;   %1/2 period
orbit.T   = 2*T12; %period

%--------------------------------------------------------------------------
% Compute the final vector of free variables for a potential continuation
% procedure.
%--------------------------------------------------------------------------
orbit.cont.gv = [v0(1) ; v0(3) ; v0(5) ; T12];

%--------------------------------------------------------------------------
% Compute the null vector of the Jacobian for a potential continuation
% procedure.
%--------------------------------------------------------------------------
orbit.cont.nv = -null(DFX);

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