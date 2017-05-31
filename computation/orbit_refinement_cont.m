function orbit = orbit_refinement_cont(cr3bp, orbit, params, cst)
% ORBIT_REFINEMENT_CONT computation of halo, vertical, and planar lyapunov 
% symmetric periodic orbits in the CRTBP, from an initial guess, as part of
% a continuation procedure.
%
% ORBIT = ORBIT_REFINEMENT_CONT(CR3BP, ORBIT, PARAMS, YVG, CST) computes an 
% orbit in the system CR3BP, with the desired characteristics
% listed in the structure ORBIT (size or energy, lagrange point), and from
% an first guess of initial condition YVG.
% A differential correction process is applied to get a real periodic 
% orbit (see diff_corr_3D_cont).
% This differential corrector is heavily based on the summary by Pavlak in
% his PhD thesis (2013), although the basic principles are much older than
% his work, and is best summarized in most of Doedel et al. papers (2003,
% 2007).
%
% See Pavlak 2013, section 3.3, for details <a href="matlab: 
% web('https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2013_Pavlak.pdf','-browser')">(link)</a>.
%
% Finally, a post process routine is applied in order to compute various
% characteristics of the orbit, listed below. Depending on the user-defined
% parameter structure PARAMS, the result can be plotted at the end of the
% routine.
%
% At the end of this routine, the following parameters are updated in the
% orbit structure:
%
%   - ORBIT.y0:  the initial conditions.
%   - ORBIT.T12: the half period.
%   - ORBIT.T:   the full period.
%   - Either the couple (Az, Azdim) - vertical extension for halo and
%   vertical orbits, or the couple (Ax, Axdim), maximum planar extension
%   for planar lyapunov orbits.
%   - ORBIT.C: the jacobian constant
%   - ORBIT.E: the energy
%   - ORBIT.yv: the state along the orbit on a given grid over the interval
%   [0, orbit.T].
%   - ORBIT.monodromy: the monodromy matrix.
%   - ORBIT.eigenvalues: the eigenvalues of the monodromy matrix in vector
%   form.
%   - ORBIT.stable_direction: the stable eigenvector of the monodromy
%   matrix.
%   - ORBIT.unstable_direction: the unstable eigenvector of the monodromy
%   matrix.
%
% See also DIFF_CORR_3D_CONT 
%
% BLB 2016

%----------------------------------------------------------------------
% Updating the state along the current tangent to the free-variables
% vector family:
%   X0^(n+1) = X0^n + ds * DX0^n
% with  - X0^n:   the converged solution of the previous step.
%       - DX0^n:  the null vector of the Jacobian associated to X0^n
%
% Here, recall that X0 = [x0 z0 vy0 T12]^T,
% with vy0 = dot(y0), T12 = half period
%----------------------------------------------------------------------
yv0    = orbit.y0;
yv0(1) = orbit.y0(1) + orbit.cont.ds*orbit.cont.nv(1);
yv0(3) = orbit.y0(3) + orbit.cont.ds*orbit.cont.nv(2);
yv0(5) = orbit.y0(5) + orbit.cont.ds*orbit.cont.nv(3);
T12    = orbit.T12   + orbit.cont.ds*orbit.cont.nv(4);

%-------------------------------------------------------------------------%
% Initial vector is concatenated with the STM matrix (Identity matrix)
% (just in case it has not been done before).
%-------------------------------------------------------------------------%
yv0 = matrixToVector(yv0, cst.orbit.STM0, 6, 6, 6);

%-------------------------------------------------------------------------%
% Differential correction with pseudo-arclength constraint.
%-------------------------------------------------------------------------%
orbit = diff_corr_3D_cont(yv0, T12, cr3bp, orbit, params, cst);

%--------------------------------------------------------------------------
% Status
%--------------------------------------------------------------------------
orbit.status = cst.orbit.REAL;

%--------------------------------------------------------------------------
% Postprocess. After this step, the following elements are updated in the
% orbit structure:
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
%--------------------------------------------------------------------------
orbit = orbit_postprocess(cr3bp, orbit, params, cst);

%--------------------------------------------------------------------------
% Plotting (potentially)
%--------------------------------------------------------------------------
if(params.plot.orbit) %plotting
    orbit_plot(orbit, params);
end

end