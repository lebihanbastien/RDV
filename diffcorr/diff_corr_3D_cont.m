function orbit = diff_corr_3D_cont(v0, T12, cr3bp , orbit, params, cst)
% DIFF_CORR_3D_CONT Differential correction to compute for 3D periodic 
% orbits of the CRTBP, symmetric with respect to the xz-plane (halo, 
% vertical). This version include a pseudo-arclength continuation module.
%
% This differential corrector is heavily based on the summary by Pavlak in
% his PhD thesis (2013), although the basic principles are much older than
% his work, and is best summarized in most of Doedel et al. papers (2003,
% 2007).
%
% See Pavlak 2013, section 3.3, for details <a href="matlab: 
% web('https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2013_Pavlak.pdf','-browser')">(link)</a>.
%
% This routines iteratively correct the initial state v0 via an iterative
% Newton method, with an additional pseudo-arclength constraint.
% The free variables X0 are:
%       X0 = [x0 z0 vy0 T12]^T, with vy0 = dot(y0), T12 = half period.
%
% The initial targeted constraint is FX = F(X) = [y vx vz]^T = 0.
% 
% Suppose that we know the solution (X0i, DX0i) from a previous iteration
% of this routine: X0i is the corrected state, and DX0i is the null vector
% of the Jacobian matrix DFX = dF(X)/dX0. Then the pseudo-arclength
% constraint is addedto the existing constraint vector FX. This constraint
% is written as:
%       (X0 - X0i)^T * DX0i - ds = 0, 
%  where ds is the arclength stepsize (user-defined in orbit.cont.ds). 
%
% Then the constraint is augmented in the form of:
%
%   GX = G(X) = |           F(X)           | = 0
%               | (X0 - X0i)^T * DX0i - ds |
%
% The Jacobian associated to this constraint vector is:
%
%   DGX = | DFX    |
%         | DX0i^T |
%
%  with still:
%
%   DFX = | Phi21   Phi23  Phi25  vy |
%	      | Phi41   Phi43  Phi45  ax |
%         | Phi61   Phi63  Phi65  az |
%
%  with vy = dot(y), ax = dot(vx) = dot(dot(x)), same for az.
%  and where Phi is the State Transition Matrix (STM), numerically
%  integrated along with the state.
%
%  Then, the correction dX0 to be applied to X0 satisfies the following equation:
%           dX0 = inv(DGX)*GX
%
% At the end of the routine, the following elements are updated:
% - orbit.y0:         initial conditions.
% - orbit.T12:        half period.
% - orbit.T:          period.
% - orbit.cont.gv:    final vector of free variables for a potential continuation procedure.
% - orbit.cont.nv:    null vector of the Jacobian for a potential continuation procedure.
% - orbit.C:          jacobi constant
% - orbit.E:          energy constant
% - orbit.cont.iter:  number of iterations in the differential corrector
% process (can monitored outside of this routine, e.g. in order to adapt
% the stepsize orbit.cont.ds.
%
% Important remark: the null vector orbit.cont.nv is forced to be in the
% same direction as its previous instance (positive dot product). This
% allows to avoid U-turn in the family, because of a change of sign in the
% routine NULL of MATLAB.
%
% BLB 2016.

%--------------------------------------------------------------------------
%Iteration counts
%--------------------------------------------------------------------------
iter = 0;

%--------------------------------------------------------------------------
% Update the MATLAB ode options
%--------------------------------------------------------------------------
if(params.computation.type == cst.computation.MATLAB)
    options = odeset('Reltol', params.ode45.RelTol, 'Abstol', params.ode45.AbsTol);
end

%--------------------------------------------------------------------------
%Diff corr loop
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
    % Update the initial vector: X0 = [x0 z0 vy0 T12]^T
    %----------------------------------------------------------------------
    X0    = [v0(1) ; v0(3) ; v0(5); T12];
    
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
        [te, ve] = ode78_cr3bp(0.0, T12, v0, 42, cr3bp.mu);
    end
    
    %----------------------------------------------------------------------
    % Update the final vector: GX = [y vx vz (X0 - X0i)^T DX0i - ds]^T.
    % The last component is the pseudo-arclength constraint.
    %----------------------------------------------------------------------
    GX     = [ve(2) ; ve(4) ; ve(6)];
    GX(4)  = (X0 - orbit.cont.gv)'*orbit.cont.nv - orbit.cont.ds;
    
    %----------------------------------------------------------------------
    % Compute the first order correction with Newton's method:
    % Targeting the constraint GX = G(X) = 0, with the free variables X0,
    % The correction dX0 to applied to these variables satisfies the
    % equation:
    %           GX = DGX*dX0
    % With DGX the matrix defined as the Jacobian:
    %           DGX = dG(X)/dX0
    % In our case, the free variables are:
    %           X0 = [x0 z0 vy0 T12]^T
    %----------------------------------------------------------------------
    % Build the Jacobian DFX = dF(X)/dX0, with X0 = [x0, z0, vy0 T12]^T
    %
    %   DGX = | Phi21   Phi23  Phi25  vy |
    %	      | Phi41   Phi43  Phi45  ax |
    %         | Phi61   Phi63  Phi65  az |
    %         |           DX0i^T         |
    %
    %   with vy = dot(y), ax = dot(vx) = dot(dot(x)), same for az.
    %   and where Phi is the State Transition Matrix (STM).
    %   finally, DX0i is the null vector of DGX(1:3,:) computed from the
    %   previous step of the continuation procedure.
    %----------------------------------------------------------------------
    % Build the Jacobian
    DGX = zeros(4,4);
    
    % Concatenate the elements of the STM
    DGX(1,1) = ve(findIndix(2,1));
    DGX(1,2) = ve(findIndix(2,3));
    DGX(1,3) = ve(findIndix(2,5));
    
    DGX(2,1) = ve(findIndix(4,1));
    DGX(2,2) = ve(findIndix(4,3));
    DGX(2,3) = ve(findIndix(4,5));
    
    DGX(3,1) = ve(findIndix(6,1));
    DGX(3,2) = ve(findIndix(6,3));
    DGX(3,3) = ve(findIndix(6,5));
    
    % Last column is [vy ax az]^T
    vep = cr3bp_derivatives_6(te, ve, cr3bp.mu);
    DGX(1,4) = vep(2);
    DGX(2,4) = vep(4);
    DGX(3,4) = vep(6);
    
    % Last line is the transposed null vector
    DGX(4,:) = orbit.cont.nv';
    
    %----------------------------------------------------------------------
    % Stops if precision is good enough, but after the Jacobian has been
    % built, so that it is updated at the end.
    %----------------------------------------------------------------------
    if(norm(GX) < params.diff_corr.precision);
        break;
    end
    
    %----------------------------------------------------------------------
    % The first order correction is dX0 = inv(DGX)*GX
    %----------------------------------------------------------------------
    dX0 = DGX\GX;
    
    %----------------------------------------------------------------------
    %Updating initial state
    %----------------------------------------------------------------------
    v0(1) = v0(1) - dX0(1);
    v0(3) = v0(3) - dX0(2);
    v0(5) = v0(5) - dX0(3);
    T12   = T12   - dX0(4);
end


%--------------------------------------------------------------------------
%Orbit update
%--------------------------------------------------------------------------
orbit.y0  = v0;    %initial condition
orbit.T12 = T12;   %1/2 period
orbit.T   = 2*T12; %period

%--------------------------------------------------------------------------
% Compute the null vector of the Jacobian for the next step of the
% continuation procedure.
% The null vector if forced to have the same direction as its previous
% instance if it exists (positive dot product).
%--------------------------------------------------------------------------
nv = null(DGX(1:3,:));
if(isfield(orbit.cont, 'nv'))
    orbit.cont.nv = sign(orbit.cont.nv'*nv)*nv;
else
    orbit.cont.nv = nv;
end

%--------------------------------------------------------------------------
% Compute the final vector of free variables for a potential continuation
% procedure.
%--------------------------------------------------------------------------
orbit.cont.gv = [v0(1) ; v0(3) ; v0(5) ; T12];

%--------------------------------------------------------------------------
% Energy
%--------------------------------------------------------------------------
orbit.C = jacobi(orbit.y0, cr3bp.mu);  %jacobi constant
orbit.E = -0.5*orbit.C;                %energy

%--------------------------------------------------------------------------
% Final iteration
%--------------------------------------------------------------------------
orbit.cont.iter = iter;

end


%--------------------------------------------------------------------------
% Return the right indix to find Phi(i0, j0) in ve
%--------------------------------------------------------------------------
function shiftedIndix = findIndix(i0, j0)
shiftedIndix = 6 + 6*(i0-1)+j0;
end