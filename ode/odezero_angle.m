function [value,isterminal,direction] = odezero_angle(t, yvar, event)
% ODEZERO_ANGLE event routine in MATLAB ODE format.
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_ANGLE(T, YVAR, EVENT) is an event
% routine in MATLAB ODE format (see event in MATLAB help). The condition of
% the event triggering is a given angle between:
%   - the line between the current state YVAR and the center EVENT.center
%   - the x-axis.
% The value of the angle that actually triggers the event is given in 
% EVENT.value.

% See also EVENT
%
% BLB 2015

%Coordinate relative to the center
xl = yvar(1) - event.center(1);
yl = yvar(2) - event.center(2);

%Current angle
phi = atan2(yl,xl);

%Event parameters
value      = phi - event.value;
isterminal = event.isterminal;
direction  = event.direction; 

end