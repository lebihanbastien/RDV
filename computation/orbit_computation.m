function orbit = orbit_computation(cr3bp, orbit, params, cst, varargin)
% ORBIT_COMPUTATION generic routine for the computation of halo, vertical
% and planar lyapunov symmetric periodic orbits in the CRTBP.
%
% ORBIT = ORBIT_COMPUTATION(CR3BP, ORBIT, PARAMS, CST) is the default call
% to compute an orbit in the system CR3BP, with the desired characteristics
% listed in the structure ORBIT (size or energy, lagrange point).
% It makes use of a third-order richardson first guess to initialize the 
% initial conditions. 
% Then, a differential correction process is applied to get a real periodic 
% orbit.
% Finally, a post process routine is applied in order to compute various
% characteristics of the orbit, listed below. Depending on the user-defined
% parameter structure PARAMS, the result can be plotted at the end of the
% routine.
%
% See Koon et al. 2006, chapter 6, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% ORBIT = ORBIT_COMPUTATION(CR3BP, ORBIT, PARAMS, CST, DIFFCORR) gives the
% orbit with a specific type differential correction procedure, defined by
% the integer constant DIFFCORR. The values allowed for DIFFCORR
% are listed in the field CST.CORR. The value DIFFCORR overwrite the
% default differential corrector defined in the PARAMS structure.
%
% ORBIT = ORBIT_COMPUTATION(CR3BP, ORBIT, PARAMS, CST, DIFFCORR, YVG)
% computes the orbit with a specific type differential correction procedure
% defined by the integer constant DIFFCORR. The values allowed for DIFFCORR
% are listed in the field CST.CORR. The value DIFFCORR overwrite the
% default differential corrector defined in the PARAMS structure.
% Moreover, a first guess YVG is provided by the user, and no third-order
% guess is used.
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
%
% BLB 2016

%-------------------------------------------------------------------------%
% Switch on the number of inputs
%-------------------------------------------------------------------------%
switch(nargin)
    case 4 % no additionnal input: varargin = {}.
        % Initialization from third order approximation
        yvg =  third_order_orbit(orbit, 0.0, cst);
        % Computation with default parameters
        orbit = orbit_refinement(cr3bp, orbit, params, yvg, cst);
    case 5 % 1 additionnal input: DIFFCORR
        yvg =  third_order_orbit(orbit, 0.0, cst);
        % has to be used: varargin = {diffcorrtype}.
        orbit = orbit_refinement(cr3bp, orbit, params, yvg, cst, varargin{1});  
    case 6 % 2 additionnal inputs: DIFFCORR and YVG
        orbit = orbit_refinement(cr3bp, orbit, params, varargin{2}, cst, varargin{1});
    otherwise
        error('Wrong number of inputs.');
end

end