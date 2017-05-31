function dx = int_T(t,x0,mu)
%
% Circular Restricted Three-body Problem (Equations of Motion):
%          For use with a variable step integrator...
%
% Input: x0 - initial state vector: x0(1:6) stato assoluto target;
%                                   x0(7:12) absolute chaser
%                                   x0(13:48) STM TARGET


% Initialize Derivative Vector:
dx = zeros(12,1);


%% target
x_T = x0(1);
y_T = x0(2);
z_T = x0(3);
xdot_T = x0(4);
ydot_T = x0(5);
zdot_T = x0(6);

x_C = x0(7);
y_C = x0(8);
z_C = x0(9);
xdot_C = x0(10);
ydot_C = x0(11);
zdot_C = x0(12);

% Phi = reshape(x0(13:48),[6,6]);


% Define distances:
d_T = sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 );  %terra
r_T = sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 );  %luna

% Integrate: dx/dt = xdot
dx(1:3) = x0(4:6);

% Integrate: dx^2/dt^2 = EOM   (absolute equations for target)
C1 = (1-mu)/d_T^3;
C2 = mu/r_T^3;

dx(4) =  2*ydot_T + x_T  -  C1*(x_T + mu) - C2*(x_T - 1 + mu);
dx(5) = -2*xdot_T + y_T  -  C1*y_T        - C2* y_T          ;
dx(6) =                  -  C1*z_T        - C2* z_T          ;


%% chaser

% Define distances:
d_C = sqrt( (x_C + mu)^2 + y_C^2 + z_C^2 );  %terra
r_C = sqrt( (x_C - 1 + mu)^2 + y_C^2 + z_C^2 );  %luna

% Integrate: dx/dt = xdot
dx(7:9) = x0(10:12);

% Integrate: dx^2/dt^2 = EOM   (absolute equations for target)
C3 = (1-mu)/d_C^3;
C4 = mu/r_C^3;

dx(10) =  2*ydot_C + x_C  -  C3*(x_C + mu) - C4*(x_C - 1 + mu);
dx(11) = -2*xdot_C + y_C  -  C3*y_C        - C4* y_C          ;
dx(12) =                  -  C3*z_C        - C4* z_C          ;


%% STM TARGET
% 
% f4x = 1 - C1 - C2 + 3*C1/(d_T^2)*(x_T+mu)^2 + 3*C2/(r_T^2)*(x_T+mu-1)^2;
% f5y = 1 - C1 - C2 + 3*y_T^2*(C1/(d_T^2)     +   C2/(r_T^2));
% f6z =   - C1 - C2 + 3*z_T^2*(C1/(d_T^2)     +   C2/(r_T^2));
% 
% f4y = 3*y_T*    (C1/(d_T^2)*(x_T+mu) +  C2/(r_T^2)*(x_T+1-mu));
% f4z = 3*z_T*    (C1/(d_T^2)*(x_T+mu) +  C2/(r_T^2)*(x_T+1-mu));
% f5z = 3*z_T*y_T*(C1/(d_T^2)          +  C2/(r_T^2));


% f4x = 1 - (1-mu)/(d_T^3) - mu/(r_T^3) + 3*(1-mu)*((x_T+mu)^2)/(d_T^5) + 3*mu*((x_T+mu-1)^2)/(r_T^5);
% f5y = 1 - (1-mu)/(d_T^3) - mu/(r_T^3) + 3*(1-mu)*     (y_T^2)/(d_T^5) + 3*mu*       (y_T^2)/(r_T^5);
% f6z =   - (1-mu)/(d_T^3) - mu/(r_T^3) + 3*(1-mu)*     (z_T^2)/(d_T^5) + 3*mu*       (z_T^2)/(r_T^5);
% 
% f4y = 3*y_T* (1-mu)*(x_T+mu)/(d_T^5) +  3*y_T*mu*(x_T+1-mu)/(r_T^5);
% f4z = 3*z_T* (1-mu)*(x_T+mu)/(d_T^5) +  3*z_T*mu*(x_T+1-mu)/(r_T^5);
% f5z = 3*y_T*z_T* (1-mu)/(d_T^5)      +  3*y_T*z_T*mu/(r_T^5);

% f1 = mu/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(3/2) - (mu - 1)/((mu + x_T)^2 + y_T^2 + z_T^2)^(3/2) - (3*mu*(2*mu + 2*x_T - 2)^2)/(4*((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2)) + (3*(2*mu + 2*x_T)^2*(mu - 1))/(4*((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2)) - 1;
% f2 = (3*y_T*(2*mu + 2*x_T)*(mu - 1))/(2*((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2)) - (3*mu*y_T*(2*mu + 2*x_T - 2))/(2*((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2));
% f3 = (3*z_T*(2*mu + 2*x_T)*(mu - 1))/(2*((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2)) - (3*mu*z_T*(2*mu + 2*x_T - 2))/(2*((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2));
% f5 =  mu/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(3/2) - (mu - 1)/((mu + x_T)^2 + y_T^2 + z_T^2)^(3/2) + (3*y_T^2*(mu - 1))/((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2) - (3*mu*y_T^2)/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2) - 1;
% f6 = (3*y_T*z_T*(mu - 1))/((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2) - (3*mu*y_T*z_T)/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2);
% f9 =  mu/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(3/2) - (mu - 1)/((mu + x_T)^2 + y_T^2 + z_T^2)^(3/2) + (3*z_T^2*(mu - 1))/((mu + x_T)^2 + y_T^2 + z_T^2)^(5/2) - (3*mu*z_T^2)/((mu + x_T - 1)^2 + y_T^2 + z_T^2)^(5/2);
% 
% A = [zeros(3)    eye(3)
%      f1 f2 f3    0 2 0
%      f2 f5 f6   -2 0 0
%      f3 f6 f9    0 0 0];
% 
% Phi_dot = A * Phi;
% 
% Phi_vec = reshape(Phi_dot,1,[]);
% 
% for i=1:length(Phi_vec)
%     dx(i+12) = Phi_vec(i);
% end

end