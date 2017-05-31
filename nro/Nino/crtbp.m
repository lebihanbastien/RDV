function dx = crtbp(t,x0,mu)
%
% Circular Restricted Three-body Problem (Equations of Motion):
%          For use with a variable step integrator...
%
% Input: x0 - initial state vector: x0(1:6) relative state in the synodic frame;
% x0(7:12) target state in the synodic frame; x0(13:48) STM in the synodic
% frame
%
%

x = x0(1);
y = x0(2);
z = x0(3);
xdot = x0(4);
ydot = x0(5);
zdot = x0(6);

x_T = x0(7);
y_T = x0(8);
z_T = x0(9);
xdot_T = x0(10);
ydot_T = x0(11);
zdot_T = x0(12);

Phi = reshape(x0(13:48),[6,6]);

% Initialize Derivative Vector:
dx = zeros(48,1);

% Define distances:
d_T = sqrt( (x_T + mu)^2 + y_T^2 + z_T^2 );  %terra
r_T = sqrt( (x_T - 1 + mu)^2 + y_T^2 + z_T^2 );  %luna

d = sqrt( (x_T + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 );  %d_T + x
r = sqrt( (x_T - 1 + mu + x)^2 + (y_T + y)^2 + (z_T + z)^2 ); %r_T + x

% Define constants
C1 = (1-mu)/d_T^3;
C2 = mu/r_T^3;
C3 = (1 - mu)/d^3;
C4 = mu/r^3;


% Relative state in SYN
% Integrate: dx/dt = xdot 
dx(1:3) = x0(4:6);

% Integrate: dx^2/dt^2 = EOM  (relative equations in synodic frame)
dx(4) =   2*ydot + x + C1*(x_T + mu) - C3*(x_T + x + mu) + C2*(x_T + mu - 1) - C4*(x_T + x + mu - 1);
dx(5) =  -2*xdot + y + C1*y_T - C3*(y + y_T) + C2*y_T - C4*(y + y_T);
dx(6) =                C1*z_T - C3*(z + z_T) + C2*z_T - C4*(z + z_T);


% Target state in SYN
% Integrate: dx/dt = xdot
dx(7:9) = x0(10:12);

% Integrate: dx^2/dt^2 = EOM   (absolute equations for target)
dx(10) =  2*ydot_T + x_T  -  C1*(x_T + mu) - C2*(x_T - 1 + mu);
dx(11) = -2*xdot_T + y_T  -  C1*y_T        - C2* y_T          ;
dx(12) =                  -  C1*z_T        - C2* z_T          ;


%State Transition Matrix
A = A_Nino([x y z], [x_T y_T z_T], mu);
Phi_dot = A * Phi;

Phi_vec = reshape(Phi_dot,1,[]);

for i=1:length(Phi_vec)
    dx(i+12) = Phi_vec(i);
end

end