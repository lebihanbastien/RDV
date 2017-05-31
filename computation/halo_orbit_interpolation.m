function orbit = halo_orbit_interpolation(cr3bp, orbit, abacus, params, cst)
% HALO_ORBIT_INTERPOLATION computes a halo orbit of the Earth-Moon system.
%
% HALO_ORBIT_INTERPOLATION(CR3BP, ORBIT, ABACUS, PARAMS, CST) computes a
% halo orbit of the Earth-Moon CRTBP, with the characteristics given in
% ORBIT:
%     - ORBIT.li: the lagrange point of reference.
%     - ORBIT.family: NORTHERN or SOUTHERN.
% and either
%     - ORBIT.C: the jacobian constant.
% or
%     - ORBIT.Az: the vertical extension.
% The routines makes use of user-defined parameters and constants defined
% in the structure PARAMS and CST, respectively. The interpolation of the
% initial conditions are made with data stored in the structure ABACUS.
% The temporary structure fit contains all the elements 
% for this interpolation.
%
% WARNING 1: the interpolation alone should give a good enough guess to
% directly produce a periodic orbit. Hence the differential correction
% process performed in the routine orbit_refinement should be very fast, or
% even useless.
%
% WARNING 2: the interpolation possibilies are limited: see energy bounds
% in halo_init.
%
% RM: the data files contain the NORTHERN family. The SOUTHERN family is
% obtained by simple symmetry with respect to the xy-plane.
%
% BLB 2015

%--------------------------------------------------------------------------
% First guess from abacus
%--------------------------------------------------------------------------
if(isfield(orbit, 'C')) %if an energy level is provided
    
    if(orbit.C > abacus.Cjaclimit(1) && orbit.C < abacus.Cjaclimit(2))
        fit.half = 2;
        fit.degree = 2*fit.half;
        [~, array_position] = min(abs(abacus.matrix(:,8) - orbit.C));
        fit.x =  abacus.matrix(array_position - fit.half: array_position + fit.half,8);
        for count =1:6
            fit.y(:,count) =  abacus.matrix(array_position - fit.half: array_position + fit.half,count);
            %Fitting for every dimension of the state (6)
            [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
            %Evaluation
            fit.f(count) = polyval(fit.p,orbit.C,[],fit.mu);
        end
    else
        disp('WARNING: the desired jacobi cst is out of bounds in halo_init.matrix');
    end
    
elseif(isfield(orbit, 'Az')) %if a vertical extension is provided
    
    if(orbit.Az < abacus.Azlimit)
        fit.half = 2;
        fit.degree = 2*fit.half;
        [~, array_position] = min(abs(abacus.matrix(:,7) - orbit.Az));
        fit.x =  abacus.matrix(array_position - fit.half: array_position + fit.half,7);
        for count =1:6
            fit.y(:,count) =  abacus.matrix(array_position - fit.half: array_position + fit.half,count);
            %Fitting for every dimension of the state (6)
            [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
            %Evaluation
            fit.f(count) = polyval(fit.p,orbit.Az,[],fit.mu);
        end
    else
        disp('WARNING: the vertical extension Az is out of range in halo_init.matrix');
    end
    
    
end

% First guess
yv0 = fit.f(1:6); 

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.family.SOUTHERN))
    yv0(3) = -yv0(3);    
end


%--------------------------------------------------------------------------
% Refinement
%--------------------------------------------------------------------------     
orbit = orbit_refinement(cr3bp, orbit, params, yv0, cst);

end