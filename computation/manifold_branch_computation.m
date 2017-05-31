function manifold_branch = manifold_branch_computation(cr3bp, orbit, manifold_branch, pos, t, params, cst, varargin)
% MANIFOLD_BRANCH_COMPUTATION computes a manifold leg associated to a given
% orbit in the CRTBP.
%
% MANIFOLD_BRANCH_COMPUTATION(CR3BP, ORBIT, MANIFOLD_BRANCH, POS, T,
% PARAMS, CST) computes a manifold leg in the system CR3BP associated to
% the orbit ORBIT, from information contained in the structure
% MANIFOLD_BRANCH. The starting position on the orbit is given by the
% variable POS (POS in [0 1] covers the entire orbit, but any value POS>0 
% is accepted). The manifold leg is then integrated up to the time T
% (resp. -T) if is unstable (resp. stable). Moreover, an event structure
% must be provided in the MANIFOLD_BRANCH structure (see 
% INIT_MANIFOLD_BRANCH). This structure defines the condition to terminate
% the integration of the manifold leg before |t| = T is reached (e.g.
% intersection with a given plane, at a given angle, etc).
%
% See Koon et al. 2006, chapter 7, for details <a href="matlab: 
% web('http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf','-browser')">(link)</a>. 
%
% MANIFOLD_BRANCH_COMPUTATION(CR3BP, ORBIT, MANIFOLD_BRANCH, POS, T,
% PARAMS, CST, EVENT) does the same thing, but the termination event is
% user-provided in the form of an event structure EVENT.
%
% MANIFOLD_BRANCH_COMPUTATION(CR3BP, ORBIT, MANIFOLD_BRANCH, POS, T,
% PARAMS, CST, EVENT_H) does the same thing, but the termination event is
% user-provided in the form of a function handler EVENT_H.
%
% See also INIT_MANIFOLD_BRANCH
%
% BLB 2016

%--------------------------------------------------------------------------
%Switch on the number of inputs
%--------------------------------------------------------------------------
switch(nargin)
    case 7
        %------------------------------------------------------------------
        %default number of inputs: the EVENT is defined from the
        %information contained in MANIFOLD_BRANCH
        %------------------------------------------------------------------
        event_type = 'event_structure';
        
        %Check that an EVENT structure was provided by the user inside the
        %MANIFOLD_BRANCH
        if(~isfield(manifold_branch, 'event'))
            error('No EVENT structure was provided inside the MANIFOLD_BRANCH.')
        end
        
    case 8 
        %------------------------------------------------------------------
        % A user-defined event has been provided.
        % Switch on the class of the user-provided event.
        %------------------------------------------------------------------
        switch(class(varargin{1}))
            
            case 'function_handle' 
                %----------------------------------------------------------
                % a user-defined specific event function
                % handler was provided.
                %----------------------------------------------------------
                event_type = 'event_handler';
                %----------------------------------------------------------
                % First, we need to check that only MATLAB routines will be
                % used, since MEX functions do not support user-defined
                % event routines.
                %----------------------------------------------------------
                if(params.computation.type ~= cst.computation.MATLAB)
                    error(['A user-defined event function was provided but MEX functions are used. ',...
                        'MEX functions do not support user-defined event routines. ',...
                        'Change the default computation type to cst.computation.MATLAB. ']);
                end
                %----------------------------------------------------------
                % If everything is ok, we get the event routine handle from
                % varargin.
                %----------------------------------------------------------
                fevent = varargin{1};
                
            case 'struct'
                %----------------------------------------------------------
                % A user-defined user structure has been provided
                % After this step, we are in a situation exactly equivalent
                % of nargin = 7. In particular, the event structure
                % MANIFOLD_BRANCH.EVENT is overwritten.
                %----------------------------------------------------------
                event_type = 'event_structure';
                manifold_branch.event = varargin{1}; 
            otherwise
                error('Wrong class for the user-provided event.');
        end
        
    otherwise
        error('Wrong number of inputs.');
end

%--------------------------------------------------------------------------
% Prepare the event if a structure has been provided
%--------------------------------------------------------------------------
if(strcmp(event_type,'event_structure'))
    % Switch on the event type in the manifold
    switch(manifold_branch.event.type)
        case cst.manifold.event.type.FREE
            %nothing to do here, no event.
            
        case cst.manifold.event.type.ANGLE_SECTION
            fevent =  @(t,y)odezero_angle(t,y,manifold_branch.event);
            
        case {cst.manifold.event.type.X_SECTION,...
                cst.manifold.event.type.Y_SECTION,...
                cst.manifold.event.type.Z_SECTION}
            fevent = @(t,y)odezero_plane(t,y,manifold_branch.event);
            
        case cst.manifold.event.type.FLIGHT_PATH_ANGLE
            fevent = @(t,y) odezero_flightpathangle(t,y,manifold_branch.event);
            
        otherwise
            error('Unknown event type.');
    end
end

%--------------------------------------------------------------------------
% Integration on the orbit until the position POS is reached.
% POS should satisfy POS >= 0.
%--------------------------------------------------------------------------
if(pos > 0)
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
        [~,yvv] = ode113(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 orbit.T*pos],orbit.y0,options);
        yv = yvv(end,:);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        [~, yv] = ode78_cr3bp(0.0, orbit.T*pos, orbit.y0, 42, cr3bp.mu);
    end
else
    yv = orbit.y0';
end

%--------------------------------------------------------------------------
% Current STM at POS
%--------------------------------------------------------------------------
STM = eye(6);
for i = 1 : 6
    for j = 1 : 6
        m = 6*(i-1) + j;
        STM(i,j) = yv(m+6);
    end
end

%--------------------------------------------------------------------------
% New vectors and initial state at POS
%--------------------------------------------------------------------------
if(manifold_branch.stability == cst.manifold.STABLE)
    vecp = STM*orbit.stable_direction;%New vector
    vecp_norm = vecp/norm(vecp);
else
    vecp = STM*orbit.unstable_direction;%New vector
    vecp_norm = vecp/norm(vecp);
end

%Initial state
xs01 = (1:6)';
for i = 1 : 6
    xs01(i) = yv(end,i) + vecp_norm(1)*manifold_branch.way*cr3bp.d_man*vecp_norm(i);
end


%--------------------------------------------------------------------------
%Integration direction
%--------------------------------------------------------------------------
% Backwards or forward integration
if(manifold_branch.stability == cst.manifold.STABLE)
    tspan = [0 -t];
else
    tspan = [0 t];
end

%--------------------------------------------------------------------------
% Termination condition
%--------------------------------------------------------------------------
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    switch(event_type)
        case 'event_structure'
            switch(manifold_branch.event.type)
                case cst.manifold.event.type.FREE
                    options = odeset('Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
                    [t,yv] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
                    te = t(end);
                    yve = yv(end, :);
                otherwise
                    options = odeset('Event',fevent, 'Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
                    [t,yv,te,yve] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
                    %If no event
                    if(isempty(te))
                        te = t(end);
                        yve = yv(end, :);
                    end
            end
        case 'event_handler'  %the user-defined routine overwrite the potential EVENT structure in MANIFOLD_BRANCH
            options = odeset('Event',fevent, 'Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);
            [t,yv,te,yve] = ode113(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu),tspan,xs01,options);
            %If no event
            if(isempty(te))
                te = t(end);
                yve = yv(end, :);
            end
    end
    
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    switch(manifold_branch.event.type)
        case cst.manifold.event.type.FREE
            [te, yve, ~, yv] = ode78_cr3bp(0.0, tspan(2), xs01, 6, cr3bp.mu);
        otherwise
            [te, yve, ~, yv] = ode78_cr3bp_event(0.0, tspan(2), xs01, 6, cr3bp.mu, manifold_branch.event);
    end
end

%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
manifold_branch.termination_time = te;
manifold_branch.yv = yve;

%--------------------------------------------------------------------------
% Plotting (potentially)
%--------------------------------------------------------------------------
if(params.plot.manifold_branch)
    manifold_plot(yv, orbit, manifold_branch, params, cst);
end


end