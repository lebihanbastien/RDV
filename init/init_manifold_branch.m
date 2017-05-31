function manifold_branch = init_manifold_branch(stability, way, varargin)
% INIT_MANIFOLD_BRANCH initializes a manifold leg.
%
% MANIFOLD_BRANCH = INIT_MANIFOLD_BRANCH(STABILITY, WAY) initializes a
% manifold leg with a given STABILITY (STABLE/UNSTABLE) and a given WAY
% (EXTERIOR/INTERIOR).
%
% MANIFOLD_BRANCH = INIT_MANIFOLD_BRANCH(STABILITY, WAY, EVENT) initializes
% a manifold leg with a given STABILITY (STABLE/UNSTABLE), a given WAY
% (EXTERIOR/INTERIOR), and a given EVENT structure (see INIT_EVENT.m).
%
% MANIFOLD_BRANCH = INIT_MANIFOLD_BRANCH(STABILITY, WAY, EVENT_TYPE,
% EVENT_VALUE, EVENT_ISTERMINAL, EVENT_DIRECTION, EVENT_CENTER, CST)
% initializes a manifold leg with a given STABILITY (STABLE/UNSTABLE),  
% a given WAY (EXTERIOR/INTERIOR), and a given EVENT structure, defined via
% the various user-provided EVENT_* parameters. The structure CST is
% required.
%
% See also INIT_EVENT
%
% BLB 2016.
        
%--------------------------------------------------------------------------
switch(nargin)
    
    case 2
        %------------------------------------------------------------------
        % if nargin == 2 we do not have any additionnal inputs in varargin.
        % There is no event defined and the manifold is let free to be
        % integrated indefinitely.
        %------------------------------------------------------------------
        %Stability (stable or unstable)
        manifold_branch.stability = stability;
        %Exterior or interior
        manifold_branch.way = way;
        %Event type during integration
        manifold_branch.event.type = 'FREE';
        
    case 3
        %------------------------------------------------------------------
        % if nargin == 3 we have the following inputs in varargin:
        % 1. event (structure).
        %------------------------------------------------------------------
        %Stability (stable or unstable)
        manifold_branch.stability = stability;
        %Exterior or interior
        manifold_branch.way = way;
        %Event during integration
        manifold_branch.event = varargin{1};
        
    case 8
        %------------------------------------------------------------------
        % if nargin == 8 we have the following inputs in varargin, in this order:
        % 1. event_type
        % 2. event_value
        % 3. event_isterminal
        % 4. event_direction
        % 5. event_center
        % 6. cst
        % The structure manifold_branch.event is built from these
        % informations.
        %------------------------------------------------------------------
        cst = varargin{6};
        %Stability (stable or unstable)
        manifold_branch.stability = stability;
        %Exterior or interior
        manifold_branch.way = way;
        %Event type during integration
        manifold_branch.event.type = varargin{1};
        %Definition parameter for the event plane
        manifold_branch.event.value = varargin{2};
        %Is the event terminal?
        manifold_branch.event.isterminal = varargin{3};
        %Dimension: dim = 1 for x, dim = 2 for y, dim = 0 if there is no event set
        switch(varargin{1})
            case cst.manifold.event.type.X_SECTION
                manifold_branch.event.dim = 1;
            case cst.manifold.event.type.Y_SECTION
                manifold_branch.event.dim = 2;
            case cst.manifold.event.type.Z_SECTION
                manifold_branch.event.dim = 3;
            case cst.manifold.event.type.FREE
                manifold_branch.event.dim = 0;
            otherwise
                manifold_branch.event.dim = 0;
        end
        
        %Decreasing (-1), increasing (1) or all (0) zeros;
        manifold_branch.event.direction = varargin{4};
        
        %Center (if the event is an angle event e.g. for lunar flyby detection).
        manifold_branch.event.center = varargin{5};
        
        %MAx events: one by default
        manifold_branch.event.max_events = 1;
    otherwise
        error('Wrong number of inputss');
end

end

