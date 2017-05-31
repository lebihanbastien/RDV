function [value,isterminal,direction] = odezero_z(t,yvar)
% ODEZERO_Z event routine in MATLAB ODE format.
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_Z(T, YVAR, EVENT) 
% is an event routine in MATLAB ODE format (see event in MATLAB help). 
% The condition of the event triggering is:
%     YVAR(3) = 0;
%
% Note that this routine is already included in ODEZERO_PLANE
%
% See also EVENT, ODEZERO_PLANE
%
% BLB 2015

%Event parameters
value = yvar(3);
isterminal = 1;
direction =  0;
end