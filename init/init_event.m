function [ event ] = init_event(type, value, isterminal, direction, center, cst)
% INIT_EVENT initializes an event structure, with the following structure:
%
%   1. event.type:  the type of event (X_SECTION, ANGLE_SECTION...)
%   2. event.value: the associated value (example: X0 for the plane X = X0)
%   3. event.isterminal: if true, the integration stops at the first event
%   (depecrated if the MEX functions are used)
%   4. event.dim:  equal to 1 for X_SECTION, 2 for Y_SECTION, etc.
%   5. event.direction: directions of crossing considered.
%   6. event.center: center for angle computation, if the event is an angle
%   section.
%   7. event.max_events: with MEX functions, any given maximum number of
%   events can be used to terminate the integration. This field is used in
%   place of event.isterminal to check wether the integration must be
%   stopped.
%
%  See also EVENT
%
%  BLB 2015

%Event type during integration
event.type = type;
%Definition parameter for the event plane
event.value = value;
%Is the event terminal?
event.isterminal = isterminal;
%Dimension: dim = 1 for x, dim = 2 for y, dim = 0 if there is no event set
switch(event.type)
    case cst.manifold.event.type.X_SECTION
        event.dim = 1;
    case cst.manifold.event.type.Y_SECTION
        event.dim = 2;
    case cst.manifold.event.type.Z_SECTION
        event.dim = 3;
    otherwise
        event.dim = 1; %by default, to avoid any error in mex files.

end

%Decreasing (-1), increasing (1) or all (0) zeros;
event.direction = direction;
%Center (if the event is an angle event e.g. for lunar flyby detection).
event.center = center;
%maximum number of events (only for mex integration)
event.max_events = 1;

end

