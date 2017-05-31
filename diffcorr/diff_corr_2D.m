function orbit = diff_corr_2D(v0, cr3bp , orbit, params, cst)
% DIFF_CORR_2D Differential correction to compute 2D planar lyapunov 
% periodic orbits of the CRTBP, symmetric with respect to the xz-plane.
%
% This routines is basically equivalent to DIFF_CORR_3D_BB but is
% restricted to the xy planar motion. The user is referred to the 3D
% routine for details.
%
% See also DIFF_CORR_3D_BB
%
% BLB 2016.

%--------------------------------------------------------------------------
%Iteration counts
%--------------------------------------------------------------------------
iter = 0;

%--------------------------------------------------------------------------
%Type of event: reaching the y = 0 plane
%--------------------------------------------------------------------------
val_par = init_event(cst.manifold.event.type.Y_SECTION,...
                     0.0,...
                     cst.manifold.event.isterminal.YES,...
                     cst.manifold.event.direction.ALL,... 
                     cr3bp.m1.pos,...
                     cst);
                 
%Options for Matlab integration
options = odeset('Events',@odezero_y,'Reltol', params.ode113.RelTol, 'Abstol', params.ode113.AbsTol);          

%--------------------------------------------------------------------------
% Differential correction loop
%--------------------------------------------------------------------------
while(true)
    iter = iter+1;
    
    if(iter > 50)
        disp('WARNING: maximum iterations reached in differential_correction');
        break;
    end
     
    %----------------------------------------------------------------------
    % Integration stops at y=0
    %----------------------------------------------------------------------
    if(params.computation.type == cst.computation.MATLAB)
        %-----------------------------
        % If MATLAB routines only
        %-----------------------------
        [~,yv,te,ve,~] = ode113(@(t,y)cr3bp_derivatives_42(t,y,cr3bp.mu),[0 10],v0,options);
    else
        %-----------------------------
        % If MEX routines are allowed
        %-----------------------------
        if(params.plot.diff_corr)
            [te, ve, ~, yv] = ode78_cr3bp_event(0.0, 10, v0, 42, cr3bp.mu, val_par);
        else
            [te, ve] = ode78_cr3bp_event(0.0, 10, v0, 42, cr3bp.mu, val_par);
        end
        
    end
    
    %----------------------------------------------------------------------
    %Update the final state: FX = [vy]
    %----------------------------------------------------------------------
    FXr = ve(4);
    
   
    %----------------------------------------------------------------------
    % Compute the first order correction with Newton's method
    %----------------------------------------------------------------------
    % If corr_type == cst.corr.Z0_FIXED:
    %                      Af                         ppf           Bf
    % [dxfpoint]  = ( [Phif41  Phif45] - 1/ypoint * [xpp] * [Phi21   Phi25] ) * [   dx0  ]
    % [dzfpoint]      [Phif61  Phif65]              [zpp]                       [dy0point]
    %
    
    % If corr_type == cst.corr.X0_FIXED:
    %                      Af                         ppf           Bf
    % [dxfpoint]  = ( [Phif43  Phif45] - 1/ypoint * [xpp] * [Phi23   Phi25] ) * [   dz0  ]
    % [dzfpoint]      [Phif63  Phif65]              [zpp]                       [dy0point]
    %
    
    %Af and Bf
    Af = ve(6+23); %Phif45
    Bf = ve(6+11); %Phif25
    
    %Derivative of ve at y=0 (t = te)
    vep = cr3bp_derivatives_6(te, ve, cr3bp.mu);
       
    %ppf vector
    ppf = vep(4)/ve(5);  %xfpp

    %New Af
    Af = Af - ppf*Bf;
     
    %Stops if precision is good enough
    if(norm(FXr) < params.diff_corr.precision);
        break;
    end
    
    %----------------------------------------------------------------------
    % The first order correction is computed as dX0 = inv(Af)*FXr
    %----------------------------------------------------------------------
    dX0 = FXr/Af;
    
    %----------------------------------------------------------------------
    %Updating initial state
    %----------------------------------------------------------------------
    v0(5) = v0(5) - dX0;  
    
    %----------------------------------------------------------------------
    % Plotting (potentially)
    %----------------------------------------------------------------------
    if(params.plot.diff_corr)
        plotDiffCorr(yv, cr3bp, iter);
    end   
   
end

%Orbit update
orbit.y0 = v0;
orbit.T12 = te; %1/2 period
orbit.T   = 2*te; %period
end

%--------------------------------------------------------------------------
% Plotting the iterations
%--------------------------------------------------------------------------
function [] = plotDiffCorr(yv, cr3bp, iter)
    figure(1)
    hold on
    if(iter==1)
        plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'g');
    else
        plot(yv(:,1)*cr3bp.L,yv(:,2)*cr3bp.L, 'r');
    end
end