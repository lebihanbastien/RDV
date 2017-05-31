function [ output, isSolution ] = lfb(manifold_branch_stable, cr3bp,  earth, moon, user, params,  cst,  varargin)
% LFB computes a LEO-to-Manifold (Halo) transfer via a differential 
% correction procedure.
%
% [OUTPUT, ISSOLUTION] = LFB(MANIFOLD_BRANCH_STABLE, CR3BP,  EARTH, MOON, 
% USER, PARAMS,  CST) takes as an input an interior stable manifold leg
% and computes a maneuver at the last state on this leg, in order to reach
% a desired LEO altitude.
%
% A first guess for this maneuver is obtained via a heuristic search. The
% procedure make sure that we actually get away from the vicinity of the
% Moon, which tends to attract the trajectories.
% Then a differential correction procedure is applied to reach the right
% LEO altitude.  The procedure is directly taken from section 4.1 of 
% Gordon's master thesis (2008). <a href="matlab: 
% web('https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Masters/2008_Gordon.pdf','-browser')">(link)</a>
%
% If a solution is found, the flag ISSOLUTION is set to true and the
% characteristics of the solutions are stored in the structure OUTPUT.
%
% [OUTPUT, ISSOLUTION] = LFB(MANIFOLD_BRANCH_STABLE, CR3BP,  EARTH, MOON, 
% USER, PARAMS,  CST, FIRSTGUESS) does the same thing, but uses a
% user-provided first guess via the array FIRSTGUESS.
% 
% BLB 2016

%--------------------------------------------------------------------------
% Init
%--------------------------------------------------------------------------
%Solution is false by default
isSolution = false;

%Arbitrary large tspan used for integration between the Moon
%neighborhood and the Earth
tspan = [0 -20];

%Maximum iterations in the Differential Correction scheme
iterMax = 50;

%--------------------------------------------------------------------------
% Completing the plot (Include the LEO)
%--------------------------------------------------------------------------

if(params.plot.XY)
    Rm1 = cr3bp.m1.Req/cr3bp.L;
    VTheta = 0:0.01:2*pi;
    X_LEO = -cr3bp.mu + (Rm1+ user.hLEOa) *cos(VTheta);
    Y_LEO =      0    + (Rm1+ user.hLEOa) *sin(VTheta);
    figure(1)
    hold on;
    plot(X_LEO*cr3bp.L, Y_LEO*cr3bp.L,':k');
end

%--------------------------------------------------------------------------
% Recovery of information at the insertion point
%--------------------------------------------------------------------------
%Last position & velocity on the manifold
output.flyby.ystate = manifold_branch_stable.yv';

%Last position & velocity on the manifold in yvcm(1:6)
output.flyby.ystatem = (1:42)';
output.flyby.ystatem(1:6) =  output.flyby.ystate;

%Inline identity matrix in yvcm(7:42)
output.flyby.ystatem = matrixToVector(output.flyby.ystatem, cst.orbit.STM0, 6, 6, 6);

%Initial velocity
vminit = output.flyby.ystate(4:6);

%Lunar Flight Path Angle (in radians)
output.flyby.lunarFlightPathAngle = flight_path_angle(output.flyby.ystate, moon.position);

%Distance to the center of the Moon
output.flyby.distanceToMoon = norm(moon.position' - output.flyby.ystate(1:3));


%--------------------------------------------------------------------------
% Correction scheme
%--------------------------------------------------------------------------

%Checking for collision with the Moon's surface prior to computation
if(output.flyby.distanceToMoon > cr3bp.m2.Rm/cr3bp.L)
    %----------------------------------------------------------------------
    % First guess for the maneuver at injection point
    % This step is crucial, and the final result depends greatly on the
    % quality of this first guess. Here, a heuristic search is performed.
    %----------------------------------------------------------------------
    if(size(varargin) == 0)
        disp('lfb. No first guess was provided by the user. A rough grid search is performed.');
        scalDV = argMinE1(earth, tspan, output, cr3bp, user, params, cst); %rough grid search
        if(user.flyby.tangentManeuver) %if tangency is forced, fminunc is used. otherwise, to much chance of failure.
            disp('lfb. Since the tangency is forced at injection point, ');
            disp('the first guess is refined with fminunc to increase the chance of success.');
            OptOptions = optimoptions('fminunc', 'TolFun', 1e-12, 'TolX', 1e-8);
            scalDV = fminunc(@(deltaV) absE1(deltaV, earth, tspan, output, cr3bp, user, params, cst), scalDV, OptOptions);
        end
        disp(['lfb. Flyby maneuver scale factor = ', num2str(scalDV)]);
        output.flyby.deltaV = output.flyby.ystatem(4:6) * scalDV;                % forced tengency
    else
        disp('lfb. A first guess was provided by the user.');
        output.flyby.deltaV = varargin{1};
    end
    
    %----------------------------------------------------------------------
    % Correction is set using the first guess
    %----------------------------------------------------------------------
    for i = 1 : 3
        output.flyby.ystatem(i+3) =  output.flyby.ystatem(i+3) + output.flyby.deltaV(i);
    end
    
    
    %----------------------------------------------------------------------
    % Heuristic preprocessing to make sure that we get away from the Moon/
    %----------------------------------------------------------------------
    iter = 0;
    while(iter < iterMax)
        
        %Concatenation of the id matrix
        output.flyby.ystatem = matrixToVector(output.flyby.ystatem, cst.orbit.STM0, 6, 6, 6);
        
        %------------------------------------------------------------------
        %Integration until a null terrestrial flight path angle is
        %reached
        %------------------------------------------------------------------
        if(params.computation.type == cst.computation.MATLAB)
            %-----------------------------
            % If MATLAB routines only
            %-----------------------------
            [~, yarc, tearc, yearc, ~] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),tspan,output.flyby.ystatem, earth.options);
        else
            %-----------------------------
            % If MEX routines are allowed
            %-----------------------------
            [tearc, yearc, ~, yarc] = ode78_cr3bp_event(0.0, -20, output.flyby.ystatem, 42, cr3bp.mu, earth.event);
        end
        
        if(isempty(yearc))  %no solution, we return
            disp('lfb. No solution has emerged from preprocessing. return.');
            return;
        end
        
        %------------------------------------------------------------------
        %Quit loop if the moon's neighborhood has been left
        %------------------------------------------------------------------
        output.flyby.distanceToMoon = norm(moon.position' - yearc(1:3)');
        if(output.flyby.distanceToMoon > 0.1)
            break;
        end
        
        %------------------------------------------------------------------
        %Otherwise, little arbitrary "push"
        %------------------------------------------------------------------
        if(user.flyby.tangentManeuver)
            output.flyby.ystatem(4:6) =  output.flyby.ystatem(4:6) * (1 + 10^(-2));  % forced tengency
        else
            output.flyby.ystatem(4) =  output.flyby.ystatem(4) * (1 + 10^(-2));      % along the x axis (arbitrary!)
        end
        
        %------------------------------------------------------------------
        %Update iter
        %------------------------------------------------------------------
        iter = iter +1;
        
        %------------------------------------------------------------------
        %Plot approximations
        %------------------------------------------------------------------
        if(user.showSteps)
            figure(1)
            hold on
            plot(yarc(:,1)*cr3bp.L ,yarc(:,2)*cr3bp.L , 'r', 'LineSmoothing','on');
        end
        
    end
    
    %----------------------------------------------------------------------
    % Differential correction process. The procedure is directly taken from
    % section 4.1 of Gordon's master thesis (2008).
    %----------------------------------------------------------------------
    iter = 0;
    while(iter < iterMax)
        %------------------------------------------------------------------
        % Computing the current error & correction
        %------------------------------------------------------------------
        %Integration until a null terrestrial flight path angle is reached
        if(params.computation.type == cst.computation.MATLAB)
            %-----------------------------
            % If MATLAB routines only
            %-----------------------------
            [~, yarc, tearc, yearc, ~] = ode45(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),tspan,output.flyby.ystatem, earth.options);
        else
            %-----------------------------
            % If MEX routines are allowed
            %-----------------------------
            [tearc, yearc, ~, yarc] = ode78_cr3bp_event(0.0, -20, output.flyby.ystatem, 42, cr3bp.mu, earth.event);
        end
        
        %------------------------------------------------------------------
        % If the event has occured, the first order correction to apply
        % to the initial state is computed.
        %------------------------------------------------------------------
        if(~isempty(yearc)) %if the event has occured
           
            %--------------------------------------------------------------
            % Update the elements at the maneuver point
            %--------------------------------------------------------------
            %Final position
            output.comp.ef = yearc(1:3)';
            %Final position wrt Earth
            output.comp.ertf = output.comp.ef + [cr3bp.mu ; 0 ; 0];
            %Final velocity
            output.comp.vf = yearc(4:6)';
            %Final velocity wrt Earth
            output.comp.vrtf = output.comp.vf;
            %Final position
            output.leo.ystate = yearc(1:6)';
            %Final STM
            STMf = vectorToMatrix(yearc, 6, 6, 6);
            
            %--------------------------------------------------------------
            % The error E1 that we seek to nullify: it is the error
            % with respect to the desired LEO altitude.
            %--------------------------------------------------------------
            E1 = norm(output.comp.ertf) - cr3bp.m1.Rm/cr3bp.L - user.hLEOa;
            
            
            %--------------------------------------------------------------
            % Check if the error is small enough to stop the procedure
            %--------------------------------------------------------------
            if(abs(E1) < params.diff_corr.precision)
                disp('lfb. A LEO solution has been found.');
                isSolution = true;
                break;
            end
            
            %--------------------------------------------------------------
            % Correction matrix:
            %
            %         | ertf(1)/|ertf|  ertf(2)/|ertf|   ertf(3)/|ertf|  |
            % Mcorr = |  vrtf(1)        vrtf(2)           vrtf(3)        |
            %         |  ertf(1)        ertf(2)           ertf(3)        |
            %
            % where ertf is the position with respect to the Earth at the
            % maneuver point, and vrtf is the corresponding velocity.
            %--------------------------------------------------------------
            Mcorr = zeros(2,6);
            Mcorr(1,1:3) = + output.comp.ertf(1:3) / norm(output.comp.ertf);
            Mcorr(2,1:3) = - output.comp.vrtf(1:3);
            Mcorr(2,4:6) = - output.comp.ertf(1:3);
            
            %--------------------------------------------------------------
            %  Concatenated matrix:
            %
            %  Merror = |Phi artf|
            %
            %  where Phi is the State Transistion Matrix (STM) at the
            %  maneuver point, and artf the acceleration at the same point.
            %--------------------------------------------------------------
            Merror = zeros(6,4);
            for i = 1 : 6
                for j = 1 : 3
                    Merror(i,j) = STMf(i,j+3);
                end
            end
            
            %Derivation of yearc
            efd = cr3bp_derivatives_6(0, yearc, cr3bp.mu);
            
            for i = 1 : 6
                Merror(i,4) = efd(i);
            end
            
            %--------------------------------------------------------------
            % Error matrix is the product Mcorr*Merror
            %--------------------------------------------------------------
            Merror_c = Mcorr * Merror;
            
            %--------------------------------------------------------------
            % Error vector
            %--------------------------------------------------------------
            vec_error = - [E1 ; 0];
            
            %--------------------------------------------------------------
            % Correction to be made on the initial conditions: the minimum
            % norm solution is selected (see  equation 4.5 in Gordon).
            %--------------------------------------------------------------
            delta_correction = Merror_c' * ( (Merror_c * Merror_c') \ vec_error);
            
            %--------------------------------------------------------------
            % If the user desires it, the tangency of the maneuver can be
            % forced (with respect to the velocity at the maneuver point).
            %--------------------------------------------------------------
            if(user.flyby.tangentManeuver)
                %Forcing the tangency of the maneuver
                [~, position_max] = max(delta_correction(1:3));
                
                switch(position_max)
                    case 1
                        delta_correction(2) = delta_correction(1) * output.flyby.ystatem(5)/output.flyby.ystatem(4);
                        delta_correction(3) = delta_correction(1) * output.flyby.ystatem(6)/output.flyby.ystatem(4);
                    case 2
                        delta_correction(1) = delta_correction(2) * output.flyby.ystatem(4)/output.flyby.ystatem(5);
                        delta_correction(3) = delta_correction(2) * output.flyby.ystatem(6)/output.flyby.ystatem(5);
                    case 3
                        delta_correction(1) = delta_correction(3) * output.flyby.ystatem(4)/output.flyby.ystatem(6);
                        delta_correction(2) = delta_correction(3) * output.flyby.ystatem(5)/output.flyby.ystatem(6);
                end
            end
            
            %--------------------------------------------------------------
            % Updating the initial conditions
            %--------------------------------------------------------------
            for i = 1 : 3
                output.flyby.ystatem(i+3) =  output.flyby.ystatem(i+3) + delta_correction(i);
            end
            
            %--------------------------------------------------------------
            %Update the iterator
            %--------------------------------------------------------------
            iter = iter + 1;
            
            %--------------------------------------------------------------
            %Plot succesive steps, if desired
            %--------------------------------------------------------------
            if(user.showSteps)
                figure(1)
                hold on
                plot(yarc(:,1)*cr3bp.L ,yarc(:,2)*cr3bp.L , 'b');
            end
            
        else  %no event, we return
            disp('lfb. Tangency at the Earth is not obtained. return.');
            isSolution = false;
            return;
        end
    end
    
    if(~isSolution)  %no solution, we return
        disp('lfb. No solution has emerged from diff corr. return.');
        return;
    end
    
    %----------------------------------------------------------------------
    %Flyby maneuver
    %----------------------------------------------------------------------
    %final velocity @ flyby maneuver
    vmfinal = output.flyby.ystatem(4:6);
    %output.flyby.deltaV (may be used to initialize the first guess for the next point
    %on the orbit
    output.flyby.deltaV = vmfinal - vminit;
    %output.flyby.deltaV [km/s]
    output.flyby.deltaV_dim = cr3bp.L/cr3bp.T*2*pi*norm(output.flyby.deltaV);
    
    
    if(abs(output.flyby.deltaV_dim) > 1) %Maneuver greater than 1 km/s
        disp('lfb. The Flyby maneuver is greater than 1 km/s. The solution is discarded');
        isSolution = false;
    end
    
    %----------------------------------------------------------------------
    %LEO maneuver
    %----------------------------------------------------------------------
    % Approximate velocity on LEO in sideral frame (centered on the Earth)
    vLEO_sid = sqrt((1-cr3bp.mu)/(cr3bp.m1.Rm/cr3bp.L + user.hLEOa));
    %Final velocity in sideral frame
    vf_sid = sqrt( norm(output.comp.vf)^2 + (output.comp.ef(1) + cr3bp.mu)^2 + output.comp.ef(2)^2 + 2*output.comp.vf(2)*(output.comp.ef(1) + cr3bp.mu) - 2 * output.comp.vf(1) * output.comp.ef(2));
    % LEO deltaV in sideral frame
    deltaVT_sid = sqrt(vLEO_sid^2 + norm(vf_sid)^2 - 2*vLEO_sid*norm(vf_sid));
    % LEO deltaV  in synodical frame (equal to deltaV_sid because the DV are equal in both frames).
    output.leo.deltaVT = deltaVT_sid;
    %output.leo.deltaVT [km/s]
    output.leo.deltaVT_dim = cr3bp.L/cr3bp.T*2*pi*output.leo.deltaVT;
    
    %----------------------------------------------------------------------
    %Total maneuver
    %----------------------------------------------------------------------
    %deltaVTot [adim]
    output.deltaV = output.leo.deltaVT + norm(output.flyby.deltaV);
    %deltaVTot [km/s]
    output.deltaV_dim = output.leo.deltaVT_dim + output.flyby.deltaV_dim;
    
    %----------------------------------------------------------------------
    %Miscellaneous parameters
    %----------------------------------------------------------------------
    %Velocity angle @ first lunar flyby
    output.flyby.lfbva = asin(norm(cross(vmfinal,vminit))/(norm(vmfinal)* norm(vminit))) * 180/pi;
    
    %Inclination with respect to Earth-Moon plan @ LEO
    output.leo.inclination = atan(output.comp.ef(3)/sqrt((output.comp.ef(1)-cr3bp.mu)^2 + output.comp.ef(2)^2))*180/pi;
    
    %Duration of the transfer
    output.Ttot = abs(manifold_branch_stable.termination_time + tearc);
    output.Ttot_dim = output.Ttot* cr3bp.T/(2*pi) * 1 / (24*3600);
    
    %----------------------------------------------------------------------
    %Final plot
    %----------------------------------------------------------------------
    strLEO   = num2str(user.hLEO);
    strDVTot = num2str(output.deltaV_dim);
    strDVM   =  num2str(output.flyby.deltaV_dim);
    strDVT   =  num2str(output.leo.deltaVT_dim);
    strTtot  = num2str(output.Ttot_dim);
    
    strTitle1 = ['Earth-to-Halo transfer from h_{LEO} = ', strLEO, ' km'];
    strTitle2 = [' Total DV: ', strDVTot, ' km/s'];
    strTitle3 = [' DV1 (LEO): ', strDVT, ' km/s'];
    strTitle4 = [' DV2 (Flyby): ', strDVM, ' km/s'];
    strTitle5 = [' Transfer time: ', strTtot, ' days'];
    
    if(params.plot.XY && isSolution)
        figure(1)
        hold on
        plot(yarc(:,1)*cr3bp.L ,yarc(:,2)*cr3bp.L , 'm', 'LineSmoothing','on');
        title({strTitle1, strTitle2, strTitle3, strTitle4, strTitle5});
    end
    
    if(params.plot.TD && isSolution)
        figure(4)
        hold on
        plot3(yarc(:,1)*cr3bp.L ,yarc(:,2)*cr3bp.L, yarc(:,3)*cr3bp.L , 'm', 'LineSmoothing','on');
        title({strTitle1, strTitle2, strTitle3, strTitle4, strTitle5});
    end
    
    
else
    if(params.plot.XY)
        figure(1)
        hold on
        title('Collision with the moon during flyby!');
    end
end




end

%--------------------------------------------------------------------------
% Subroutines
%--------------------------------------------------------------------------
function [absE1] = absE1(deltaV, earth, tspan, output, cr3bp, user, params, cst)
% ABSE1 Gives the absolute norm of the error E1 (error wrt desired LEO 
% altitude) for a flyby maneuver of the following type :
% DV = deltaV * output.manifold.ystatem(4:6)

%Correction is set
for i = 4 : 6
    output.flyby.ystatem(i) =  output.flyby.ystatem(i) + output.flyby.ystatem(i)*deltaV;
end

%Integration until a null terrestrial flight path angle is
%reached
if(params.computation.type == cst.computation.MATLAB)
    %-----------------------------
    % If MATLAB routines only
    %-----------------------------
    [~, ~, ~, yearc, ~] = ode45(@(t,y)cr3bp_derivatives_6(t,y,cr3bp.mu), tspan, output.flyby.ystatem(1:6), earth.options);
else
    %-----------------------------
    % If MEX routines are allowed
    %-----------------------------
    [~, yearc, ~, ~] = ode78_cr3bp_event(0.0, tspan(2), output.flyby.ystatem(1:6), 6, cr3bp.mu, earth.event);
end

%If a solution has been found, compute the error
if(~isempty(yearc))
    %Final position
    output.comp.ef = yearc(1:3)';
    %Final position wrt Earth
    output.comp.ertf = output.comp.ef + [cr3bp.mu ; 0 ; 0];
    
    %Error E1
    E1 = norm(output.comp.ertf) - cr3bp.m1.Rm/cr3bp.L - user.hLEOa;
    
    %Norm of error
    absE1 = abs(E1);
else
    absE1 = inf; %infinity is output
end

end


function argMinE1 = argMinE1(earth, tspan, output, cr3bp, user, params, cst)
% Gives the minimum argument deltaV on a rough grid for min{E1(deltaV)}
argMinE1 = 0.0;
MinE1 = absE1(argMinE1, earth, tspan, output, cr3bp, user, params, cst);
for di = 0.1:0.1:0.5
    E1 = absE1(di, earth, tspan, output, cr3bp, user, params, cst);
    if(E1 < MinE1)
        argMinE1 = di;
        MinE1 = E1;
    end
end
end