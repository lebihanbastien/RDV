function orbit = orbit_refinement(cr3bp, orbit, params, yvg, cst, varargin)
% ORBIT_REFINEMENT computation of halo, vertical, and planar lyapunov 
% symmetric periodic orbits in the CRTBP, from an initial guess.
%
% ORBIT = ORBIT_REFINEMENT(CR3BP, ORBIT, PARAMS, YVG, CST) computes an 
% orbit in the system CR3BP, with the desired characteristics
% listed in the structure ORBIT (size or energy, lagrange point), and from
% an first guess of initial condition YVG.
% A differential correction process is applied to get a real periodic 
% orbit.
% Finally, a post process routine is applied in order to compute various
% characteristics of the orbit, listed below. Depending on the user-defined
% parameter structure PARAMS, the result can be plotted at the end of the
% routine.
%
% See Koon et al. 2006, chapter 6, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>.
%
% ORBIT = ORBIT_REFINEMENT(CR3BP, ORBIT, PARAMS, YVG, CST, DIFFCORR)
% computes the orbit with a specific type differential correction procedure
% defined by the integer constant DIFFCORR. The values allowed for DIFFCORR
% are listed in the field CST.CORR. The value DIFFCORR overwrite the
% default differential corrector defined in the PARAMS structure.
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
% BLB 2016

%--------------------------------------------------------------------------
% Switch on the number of inputs:
% The type of the differential corrector used by this routine can be
% user-defined in the variable-size input varargin. If varargin = {} (no
% additionnal input is provided), the routine makes use of the default
% corrector type in the params structure.
%--------------------------------------------------------------------------
switch(nargin)
    case 5 %no type of differential corrector is provided, the default one is used (in the params structrure).
        diffcorr = params.diff_corr.type;
    case 6 %a differential corrector was provided by the user.
        diffcorr = varargin{1};
    otherwise
        error('Wrong number of inputs.');
end

%--------------------------------------------------------------------------
% Initialisation of the initial conditions.
%--------------------------------------------------------------------------
% Integration vector
yv0 = (1:42)';
% 6-dim state
yv0(1:6) = yvg(1:6);
% STM concatenation after the 6-dim state
yv0 = matrixToVector(yv0, cst.orbit.STM0, 6, 6, 6);

%--------------------------------------------------------------------------
% Differential correction procedure (DCP). At the end of this procedure,
% the following elements (at least) are updated in the orbit structure:
%   - orbit.y0:         initial conditions.
%   - orbit.T12:        half period.
%   - orbit.T: yv0         period.
%--------------------------------------------------------------------------
%Switch between orbit types (Halo, Vertical, Lyapunov)
switch(orbit.type)
    case cst.orbit.type.HALO
        %---------------------------------------------
        % Select the type of differential corrector procedure (DCP)
        %---------------------------------------------
        switch(diffcorr)
            case cst.corr.MIN_NORM
                %---------------------------------------------
                % In this case a 4-dimensionnal DCP is selected. The free
                % variables are X0 = [x0 z0 vy0 T12]^T
                %---------------------------------------------
                orbit = diff_corr_3D_full(yv0, cr3bp , orbit, params, cst);
                
            case {cst.corr.X0_FIXED, cst.corr.Z0_FIXED}
                %----------------------------------------------------------
                % Warning if the DCP type is wrong for big orbits.
                % If no DCP type was provided by the user (varargin = {}),
                % then the right corrector type is forced. If the DCP type
                % was directly user-provided varargin, we consider that the
                % user knows what he's doing, and a simple warning is sent.
                % The definition of 'big orbits' is loosesly based on
                % results from the Earth-Moon system.
                %----------------------------------------------------------
                if(isfield(orbit, 'Az_estimate') && orbit.Az_estimate > 0.0520 && diffcorr ~= cst.corr.X0_FIXED)
                    switch(nargin)
                        case 5 %no type of differential corrector is provided, the default one is used (in the params structrure).
                            diffcorr = cst.corr.X0_FIXED;
                        case 6 %a differential corrector was provided by the user.
                            warning(['You are trying to compute a big halo orbit. '...
                                'For these orbits, the differential corrector of type '...
                                'cst.corr.X0_FIXED is preferable.']);
                    end
                end
                
                %----------------------------------------------------------
                % In this case, a simple 3-dimensional DCP is selected. One
                % dimension is fixed (either x0 or z0), and the free
                % variables are either x0 or z0, and vy0.
                %----------------------------------------------------------
                % 1. The right dimensions are selected
                if(isequal(diffcorr, cst.corr.X0_FIXED))
                    xi0 = [3, 5];  %(z0, vy0) are corrected
                elseif(isequal(diffcorr, cst.corr.Z0_FIXED))
                    xi0 = [1, 5];  %(x0, vy0) are corrected
                end
                xif = [4, 6]; %(vx, vz) = 0 is targeted
                
                % 2. Perform the DCP
                orbit = diff_corr_3D_bb(yv0, cr3bp , orbit, xi0, xif, params, cst);
            otherwise
                error('Unknown differential corrector type.');
        end
        
    case cst.orbit.type.VLYAP
        %------------------------------------------------------------------
        % Select the type of differential corrector procedure (DCP)
        %------------------------------------------------------------------
        switch(diffcorr)
            case cst.corr.MIN_NORM
                %----------------------------------------------------------
                % In this case a 4-dimensionnal DCP is selected. The free
                % variables are X0 = [x0 z0 vy0 T12]^T
                %----------------------------------------------------------
                orbit = diff_corr_3D_full(yv0, cr3bp , orbit, params, cst);
                
            case {cst.corr.X0_FIXED, cst.corr.Z0_FIXED}
                %----------------------------------------------------------
                % In this case, a simple 3-dimensional DCP is selected. One
                % dimension is fixed (either x0 or z0), and the free
                % variables are either x0 or z0, and vy0.
                %----------------------------------------------------------
                % 1. The right dimensions are selected
                if(isequal(diffcorr, cst.corr.X0_FIXED))
                    xi0 = [3, 5];  %(z0, vy0) are corrected
                elseif(isequal(diffcorr, cst.corr.Z0_FIXED))
                    xi0 = [1, 5];  %(x0, vy0) are corrected
                end
                xif = [3, 4]; %(z, vx) = 0 is targeted
                
                % 2. Perform the DCP
                orbit = diff_corr_3D_bb(yv0, cr3bp , orbit, xi0, xif, params, cst);
            otherwise
                error('Unknown differential corrector type.');
        end
        
    case cst.orbit.type.PLYAP
        %------------------------------------------------------------------
        % In this case, a simple 2-dimensional DCP is selected. The free
        % variables are x0 and vy0.
        %------------------------------------------------------------------
        orbit = diff_corr_2D(yv0, cr3bp, orbit, params, cst);
end



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