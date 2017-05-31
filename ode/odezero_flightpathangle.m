function [value,isterminal,direction] = odezero_flightpathangle(t, yvar, event)
% ODEZERO_FLIGHTPATHANGLE event routine in MATLAB ODE format.
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_FLIGHTPATHANGLE(T, YVAR, EVENT) 
% is an event routine in MATLAB ODE format (see event in MATLAB help). 
% The condition of the event triggering is a given flight path angle of the
% current velocity in YVAR with respect to the center EVENT.center. The
% value of the flight path angle that actually triggers the event is given
% in EVENT.value.
%
% See also EVENT
%
% BLB 2015

% Position with respect to the center of event
xrt = yvar(1:3) - event.center';

% Velocity wrt the event
vrt = yvar(4:6)';

%Event parameters
value      = - vrt * xrt - event.value;
isterminal = event.isterminal;
direction  = event.direction;

end
