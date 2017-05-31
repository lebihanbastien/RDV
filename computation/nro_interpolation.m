function orbit = nro_interpolation(cr3bp, orbit, abacus, params, cst, varargin)
% NRO_INTERPOLATION computes Near Rectilinear Orbits (NRO)in the Earth-Moon
% Circular Restricted Three-Body Problem.
%
% ORBIT = NRO_INTERPOLATION(CR3BP, ORBIT, ABACUS, PARAMS, CST, NAME, VALUE)
% computes the NRO whose characteristics are defined inside the ORBIT
% structure:
%       - ORBIT.li: Lagrange point of reference
%       - ORBIT.family: NORTHERN or SOUTHERN family
% The routines makes use of user-defined parameters and constants defined
% in the structure PARAMS and CST, respectively. The interpolation of the
% initial conditions are made with data stored in the structure ABACUS.
% The temporary structure fit contains all the elements 
% for this interpolation.
%
% Different ways to uniquely define the orbit are available, via the couple
% (NAME, VALUE). Namely:
%       - NAME = 'altitudeOfPerigee', allows to defined the orbit via the
%       altitude of its perigee wrt the Moon surface (in km).
%       - NAME = 'perigeeDistanceToMoonCenter', allows to defined the
%       via the distance of its perigee wrt to the Moon center. We
%       basically have:
%     altitudeOfPerigee = perigeeDistanceToMoonCenter - Moon's radius.
%       - NAME = 'Az', allows to defined the orbit via its vertical
%       extension, as it is usually done with Halo orbits.
%
% WARNING 1: the interpolation alone should give a good enough guess to
% directly produce a periodic orbit. Hence the differential correction
% process performed in the routine orbit_refinement should be very fast, or
% even useless.
%
% WARNING 2: the interpolation possibilies are limited: see energy bounds
% in the abacus.
%
% RM: the data files contain the NORTHERN family. The SOUTHERN family is
% obtained by simple symmetry with respect to the xy-plane.
%
% BLB 2016

%--------------------------------------------------------------------------
% Preliminary checks
%--------------------------------------------------------------------------
if(nargin ~= 7)
    error('Wrong number of arguments');
end

if(~strcmp(cr3bp.name,'EARTH+MOON'))
    error('Wrong CRTBP. Only EARTH+MOON system is authorized.');
end

%--------------------------------------------------------------------------
% First guess from abacus
%--------------------------------------------------------------------------
switch(varargin{1})
    
    case 'altitudeOfPerigee'
        
        altitudeOfPerigee = varargin{2}/cr3bp.L;
        if(altitudeOfPerigee > abacus.altitudeOfPerigeeLimit(1) &&...
                altitudeOfPerigee < abacus.altitudeOfPerigeeLimit(2))
            %Interpolate
            fit = interpol(abacus, abacus.altitudeOfPerigee, altitudeOfPerigee);
        else
            error('The desired altitude of perigee is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.altitudeOfPerigeeLimit(1)*cr3bp.L, abacus.altitudeOfPerigeeLimit(2)*cr3bp.L);
        end
        
    case 'perigeeDistanceToMoonCenter'
        
        perigeeDistanceToMoonCenter = varargin{2}/cr3bp.L;
        if(perigeeDistanceToMoonCenter > abacus.perigeeDistanceToMoonCenterLimit(1)...
                && perigeeDistanceToMoonCenter < abacus.perigeeDistanceToMoonCenterLimit(2))
            %Interpolate
            fit = interpol(abacus, abacus.perigeeDistanceToMoonCenter, perigeeDistanceToMoonCenter);
        else
            error('The desired distance of perigee to Moon center is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.perigeeDistanceToMoonCenterLimit(1)*cr3bp.L, abacus.perigeeDistanceToMoonCenterLimit(2)*cr3bp.L);
        end
    case 'Az'
        Az = varargin{2}/cr3bp.L;
        if(Az > abacus.AzLimit(1) && Az < abacus.AzLimit(2))
            %Interpolate
            fit = interpol(abacus, abacus.Az, Az);
        else
            error('The desired distance of perigee to Moon center is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.AzLimit(1)*cr3bp.L, abacus.AzLimit(2)*cr3bp.L);
        end
        
end

% First guess
yv0 = fit.f(1:6);

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.family.SOUTHERN))
    yv0(3) = -yv0(3);
end


%--------------------------------------------------------------------------
% Refinement (differential correction)
%--------------------------------------------------------------------------
orbit = orbit_refinement(cr3bp, orbit, params, yv0, cst);

end

%--------------------------------------------------------------------------
% Interpolate the initial condition y0 to get data(y0) = value.
%--------------------------------------------------------------------------
function fit = interpol(abacus, data, value)

fit.half = 2;
[~, array_position] = min(abs(data - value));
mini = max(array_position - fit.half, 1);
maxi = min(array_position + fit.half, length(data));
fit.degree = length(mini:maxi)-1;
fit.x =  data(mini:maxi)';
for count =1:6
    fit.y(:,count) =  abacus.initialConditions(mini:maxi,count);
    %Fitting for every dimension of the state (6)
    [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
    %Evaluation
    fit.f(count) = polyval(fit.p,value,[],fit.mu);
end

end