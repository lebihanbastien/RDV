function [value,isterminal,direction] = odezero_plane(t, yvar, event)
% ODEZERO_PLANE event routine in MATLAB ODE format.
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_PLANE(T, YVAR, EVENT) 
% is an event routine in MATLAB ODE format (see event in MATLAB help). 
% The condition of the event triggering is a given value of one dimension 
% in YVAR. This dimension is given by EVENT.DIM. 
% The value of YVAR(EVENT.DIM) that actually triggers the event is given by
% the scalar EVENT.value.
%
% See also EVENT
%
% BLB 2015

%Event parameters
value = yvar(event.dim)-event.value;
isterminal = event.isterminal;
direction =  event.direction;

end